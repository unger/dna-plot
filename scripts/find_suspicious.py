"""Suggest suspicious BOLD records for review: data/bold/suspicious_candidates.tsv.

Usage: python scripts/find_suspicious.py   (reads data/bold/bold.sqlite, changes nothing)
Writes suspicious_candidates.tsv (new candidates) and suspicious_recheck.tsv (the already marked records
judged again with the current alignment and checks: 'unmark?' where no check supports the mark any more).

Checks (COI-5P records not already marked suspicious):
  bin_outlier   a record far from every other record in its BIN (> BIN_OUTLIER %), after repairing a
                single indel; BOLD's BINs are normally within ~2.2 %
  name_outlier  a named record that is clearly nearer to records of another species than to any record of
                its own species (by NAME_MARGIN %), i.e. probably misidentified or mixed up
  bin_minority_name  a named record whose BIN is dominated by another species (its own name is a small
                minority there) and whose sequence is clearly nearer that species than its own name's records
                in other BINs. Judged per record, so sequences shared by several species and species with
                only 2-3 records are not hidden
  name_swap     an X-named record in a BIN of Y plus a Y-named record in a BIN of X: the names look exchanged
  stop_codon    translation (invertebrate mitochondrial code) of the record's frame has stop codons
Review the list, then mark the ones you agree with:
  UPDATE records SET suspicious=1, suspicious_reason='...' WHERE record_id='...';
and rebuild the reference files.
"""
import collections
import itertools
import os
import re
import sqlite3

import numpy as np

import build_bin_reference as bb
import build_reference as br

BIN_OUTLIER = 3.0   # % to the nearest other record in the BIN
MISSPELL_MIN = 3     # a name must have this many records to be what another name is a misspelling of
NAME_MARGIN = 1.5   # % the nearest other-species record must be closer than the nearest same-species one
MAJORITY_SHARE = 0.6    # a BIN's majority name is at least this share of its named records
MINORITY_SHARE = 0.2    # bin_minority_name: the record's name is at most this share
SWAP_SHARE = 0.5        # name_swap: looser, for small BINs (e.g. 1 of 3)
ESTABLISHED = 3         # name_outlier: a name with this many records in a BIN, at least half of its named ones, is not a slip
MIN_COMPARED = 400
OUT = "data/bold/suspicious_candidates.tsv"
RECHECK = "data/bold/suspicious_recheck.tsv"  # the marked records, judged again
SYNONYMS = "scripts/name_synonyms.tsv"  # old name <TAB> accepted name; compared as the accepted name
STOPS = {"TAA", "TAG"}  # TGA codes for Trp in invertebrate mitochondria


def pdist(a, B):
    """% difference of code array a to each row of B over positions where both are A/C/G/T (nan if too few)."""
    ok = (a < 4) & (B < 4)
    n = ok.sum(axis=1)
    d = ((a != B) & ok).sum(axis=1)
    return np.where(n >= MIN_COMPARED, 100.0 * d / np.maximum(n, 1), np.nan)


def repaired_dist(a, b):
    r = br.repair_indels(a, b)
    ok = (r < 4) & (b < 4)
    n = ok.sum()
    return 100.0 * ((r != b) & ok).sum() / n if n >= MIN_COMPARED else np.nan


def proper_name(n):
    """Only "Genus species" counts as a name; "sp.", "group", codes and the like are placeholders."""
    return bool(re.fullmatch(r"[A-Z][a-z]+ [a-z]+", n or "")) and not n.endswith(" sp")


def edit_distance(a, b):
    prev = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        cur = [i]
        for j, y in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[-1] + 1, prev[j - 1] + (x != y)))
        prev = cur
    return prev[-1]


def name_fixes(rows):
    """old name -> (accepted name, reason). Names in SYNONYMS count as their accepted name, and a name that is
    1-2 letters from a much commoner name in the same BIN is taken as a misspelling of it. Such records are
    not suspicious, so they are compared under the accepted name."""
    fixes = {}
    if os.path.exists(SYNONYMS):
        for line in open(SYNONYMS, encoding="utf8"):
            if line.strip() and not line.startswith("#"):
                old, new = line.rstrip("\n").split("\t")[:2]
                fixes[old] = (new, "synonym")
    count = collections.Counter(r[1] for r in rows if proper_name(r[1]))
    in_bin = collections.defaultdict(set)
    for r in rows:
        if r[2] and proper_name(r[1]):
            in_bin[r[2]].add(r[1])
    for names in in_bin.values():
        for a, b in itertools.permutations(sorted(names), 2):
            if a not in fixes and count[a] * 3 <= count[b] and edit_distance(a, b) <= 2:
                fixes[a] = (b, "misspelling")
    # A misspelled name often sits in a BIN of another species (that is the interesting case), so it is
    # also matched against every name, not just those of its own BIN: 1 letter from a name with at least
    # MISSPELL_MIN records, and no other name that close.
    names = [n for n in count if n not in fixes]
    for a in names:
        close = [b for b in count if b != a and abs(len(a) - len(b)) <= 1 and edit_distance(a, b) <= 1]
        if len(close) == 1 and count[close[0]] >= MISSPELL_MIN and count[a] * 3 <= count[close[0]]:
            fixes[a] = (close[0], "misspelling")
    # follow chains (a misspelling of a synonym)
    for a in list(fixes):
        seen = {a}
        while fixes[a][0] in fixes and fixes[a][0] not in seen:
            seen.add(fixes[a][0])
            fixes[a] = (fixes[fixes[a][0]][0], fixes[a][1])
    return fixes


def stop_candidates(rows):
    """record_id -> number of stop codons, for records with stops in every reading frame. Records that
    BOLD did not pad start at varying codon positions, so all three frames are tried."""
    out = {}
    for rid, sp, b, co, nuc in rows:
        best = min(sum(1 for i in range(off, len(nuc) - 2, 3) if nuc[i:i + 3] in STOPS) for off in range(3))
        if best:
            out[rid] = best
    return out


def analyse(con, include_marked, quiet=False, marked_as_evidence=True):
    """Run all checks. include_marked=False: only records not marked suspicious (the review list).
    include_marked=True: the marked ones are judged too, as if they were new records, but they never count
    as evidence for others or in a BIN's shares (see recheck_marked), unless marked_as_evidence is set: then
    they are ordinary records. Returns (candidates, meta)."""
    marked = set()
    if include_marked and not marked_as_evidence:
        marked = {r[0] for r in con.execute("SELECT record_id FROM records WHERE suspicious=1")}
    rows = con.execute(
        """SELECT r.record_id, r.species, r.bin_uri, r.country, s.nuc FROM records r
           JOIN sequences s USING(record_id)
           WHERE r.marker_code='COI-5P' AND s.nuc IS NOT NULL""" + ("" if include_marked else " AND r.suspicious=0")).fetchall()
    rows = [(rid, sp or "", b or "", co or "", nuc.upper()) for rid, sp, b, co, nuc in rows if len(nuc) <= br.MAX_LEN]
    fixes = name_fixes(rows)
    for old, (new, why) in sorted(fixes.items()):
        if not quiet:
            print("%s -> %s (%s, %d record(s))" % (old, new, why, sum(1 for r in rows if r[1] == old)))
    rows = [(rid, fixes.get(sp, (sp,))[0], b, co, nuc) for rid, sp, b, co, nuc in rows]
    cands = collections.defaultdict(list)  # record_id -> [(check, detail)]
    meta = {r[0]: r for r in rows}

    # --- stop codons
    for rid, n in stop_candidates(rows).items():
        cands[rid].append(("stop_codon", "%d stop codon(s) in the best reading frame" % n))

    # --- put every record in its BIN's frame (BOLD does not pad all records), so columns can be compared
    glob = bb.pick_anchor([r[4] for r in rows], None)
    if glob is None:
        glob = br.encode(collections.Counter(r[4] for r in rows).most_common(1)[0][0])

    def anchor_of(seqs):
        return bb.pick_anchor(seqs, glob)
    by_bin_rows = collections.defaultdict(list)
    for r in rows:
        by_bin_rows[r[2]].append(r)
    anchors = {b: anchor_of([r[4] for r in rs]) if b else glob for b, rs in by_bin_rows.items()}
    rows = [(rid, sp, b, co, bb.in_bin_frame(nuc, anchors[b])) for rid, sp, b, co, nuc in rows]

    # --- unique sequences
    useqs = list(dict.fromkeys(r[4] for r in rows))
    idx = {s: i for i, s in enumerate(useqs)}
    R = np.full((len(useqs), br.MAX_LEN), br.PAD, dtype=np.uint8)
    for i, s in enumerate(useqs):
        R[i, :len(s)] = br.encode(s)
    recs_of = collections.defaultdict(list)  # unique sequence index -> record rows
    for r in rows:
        recs_of[idx[r[4]]].append(r)
    # only proper species names count ("Genus species"); "sp.", "group", codes and the like are placeholders
    names_of = [{r[1] for r in recs_of[i] if proper_name(r[1])} for i in range(len(useqs))]
    ref_names_of = [{r[1] for r in recs_of[i] if proper_name(r[1]) and r[0] not in marked} for i in range(len(useqs))]

    # --- BIN outliers: nearest other record in the same BIN
    by_bin = collections.defaultdict(set)
    ref_by_bin = collections.defaultdict(set)
    for r in rows:
        if r[2]:
            by_bin[r[2]].add(idx[r[4]])
            if r[0] not in marked:
                ref_by_bin[r[2]].add(idx[r[4]])
    for b, members in by_bin.items():
        members = sorted(members)
        ref_members = sorted(ref_by_bin[b])
        if len(ref_members) < 4:
            continue
        for i in members:
            others = [j for j in ref_members if j != i]
            d = pdist(R[i], R[others])
            order = np.argsort(np.where(np.isnan(d), 1e9, d))[:5]
            best = np.nanmin(d) if not np.all(np.isnan(d)) else np.nan
            if np.isnan(best) or best <= BIN_OUTLIER:
                continue
            best = min(repaired_dist(R[i], R[others[k]]) if not np.isnan(d[k]) else 1e9 for k in order)
            if best > BIN_OUTLIER:
                for r in recs_of[i]:
                    cands[r[0]].append(("bin_outlier", "%.1f %% from the nearest of %d others in %s" % (best, len(others), b)))

    # --- name outliers: nearest same-species vs nearest other-species record
    n_name_bin = collections.Counter((r[2], r[1]) for r in rows if r[2] and proper_name(r[1]) and r[0] not in marked)
    n_bin = collections.Counter(r[2] for r in rows if r[2] and proper_name(r[1]) and r[0] not in marked)

    def established(r):
        """Several records with the same name in the BIN, at least half of its named records: more than a slip,
        so the other-named records there are the odd ones (or a species COI barely separates)."""
        return n_name_bin[r[2], r[1]] >= ESTABLISHED and 2 * n_name_bin[r[2], r[1]] >= n_bin[r[2]]

    mask = collections.defaultdict(lambda: np.zeros(len(useqs), dtype=bool))
    for j, nm in enumerate(ref_names_of):
        for n in nm:
            mask[n][j] = True
    named = np.array([bool(nm) for nm in ref_names_of])
    for start in range(0, len(useqs), 1):
        for i in range(start, start + 1):
            if not names_of[i]:
                continue
            d = pdist(R[i], R)
            d[i] = np.nan
            same = np.zeros(len(useqs), dtype=bool)
            for n in names_of[i]:
                same |= mask[n]
            other = named & ~same
            ds = np.where(same & ~np.isnan(d), d, np.inf)
            do = np.where(other & ~np.isnan(d), d, np.inf)
            if not np.isfinite(do).any() or do.min() == np.inf:
                continue
            bs, bo = ds.min(), do.min()
            if bs - bo < NAME_MARGIN:
                continue
            # confirm with indel repair against the few nearest of each side
            js = [j for j in np.argsort(ds)[:5] if np.isfinite(ds[j])]
            jo = [j for j in np.argsort(do)[:5] if np.isfinite(do[j])]
            rs = min([repaired_dist(R[i], R[j]) for j in js] or [np.inf])
            ro = min([repaired_dist(R[i], R[j]) for j in jo] or [np.inf])
            rs = np.inf if np.isnan(rs) else rs
            ro = np.inf if np.isnan(ro) else ro
            # a name without any other record is no evidence (synonym, or a species COI cannot tell apart)
            if np.isfinite(rs) and rs - ro >= NAME_MARGIN:
                near = ", ".join(sorted(ref_names_of[jo[int(np.nanargmin([repaired_dist(R[i], R[j]) for j in jo]))]]))
                own = "%.1f %% from the nearest of its own name" % rs
                for r in recs_of[i]:
                    if established(r):
                        continue
                    cands[r[0]].append(("name_outlier", "%s; %.1f %% from %s" % (own, ro, near)))
    # --- a name in a BIN of another species, and exchanged names. Judged per record rather than per unique
    # sequence, so a sequence shared by several species does not hide a record's own name.
    seqs_of = collections.defaultdict(lambda: collections.defaultdict(set))  # name -> bin -> unique sequences
    named_in = collections.defaultdict(list)  # bin -> records with a proper name
    ref_named_in = collections.defaultdict(list)  # the same without the marked records
    for r in rows:
        if proper_name(r[1]):
            if r[0] not in marked:
                seqs_of[r[1]][r[2]].add(idx[r[4]])
            if r[2]:
                named_in[r[2]].append(r)
                if r[0] not in marked:
                    ref_named_in[r[2]].append(r)
    near_cache = {}

    def nearest(i, js, key):
        """% from sequence i to the nearest of the sequences js (the lower of plain and indel-repaired)."""
        if (i, key) not in near_cache:
            js = sorted(js)
            out = np.inf
            if i in js:
                out = 0.0
            elif js:
                d = pdist(R[i], R[js])
                if not np.all(np.isnan(d)):
                    order = [k for k in np.argsort(np.where(np.isnan(d), 1e9, d))[:5] if not np.isnan(d[k])]
                    out = min([np.nanmin(d)] + [x for x in (repaired_dist(R[i], R[js[k]]) for k in order) if not np.isnan(x)])
            near_cache[i, key] = out
        return near_cache[i, key]

    dissent = []  # (record, own name, majority name, bin, majority count, named count, own count, dM, dX)
    for b, mem_all in named_in.items():
        mem = ref_named_in[b]  # shares are counted over the unmarked records
        if not mem:
            continue
        cnt = collections.Counter(r[1] for r in mem)
        top = cnt.most_common(2)
        M, mc = top[0]
        if len(mem) < 3 or (len(top) > 1 and top[1][1] == mc) or mc / len(mem) < MAJORITY_SHARE:
            continue
        for r in mem_all:
            X = r[1]
            mine = 1 if r[0] in marked else 0  # a marked record is not in the counts yet
            if X == M or cnt[X] >= mc or (cnt[X] + mine) / (len(mem) + mine) > SWAP_SHARE:
                continue
            i = idx[r[4]]
            elsewhere = set().union(*[v for k, v in seqs_of[X].items() if k != b])
            dM = nearest(i, seqs_of[M][b], ("in", M, b))
            dX = nearest(i, elsewhere, ("out", X, b))
            if dM <= BIN_OUTLIER and np.isfinite(dX) and dX - dM >= NAME_MARGIN:
                dissent.append((r, X, M, b, mc, len(mem) + mine, cnt[X] + mine, dM, dX))
    by_names = collections.defaultdict(list)
    for t in dissent:
        by_names[t[1], t[2]].append(t)
    for r, X, M, b, mc, n, nx, dM, dX in dissent:
        partners = sorted(p[0][0] for p in by_names.get((M, X), []) if p[3] != b)
        if partners:
            cands[r[0]].append(("name_swap", "name probably exchanged with %s (%s sits in a %s BIN)" % (", ".join(partners[:3]), M, X)))
        if nx / n <= MINORITY_SHARE:
            cands[r[0]].append(("bin_minority_name", "%s is %d of %d named in %s (mostly %s); %.1f %% from %s, %s" % (
                X, nx, n, b, M, dM, M, "%.1f %% from its own name elsewhere" % dX)))
    return cands, meta


def recheck_marked(con, cands_all, meta_all):
    """Records marked suspicious, judged again by the current checks (the marked ones included in the
    reference, so they can support each other). A mark that no check supports any more is probably a
    leftover from an older, worse alignment: listed in RECHECK as 'unmark?'. Marks still supported are
    'keep'. The checks cannot judge records with few relatives, those come out as unmark? too, so read
    the detail before unmarking."""
    marked = con.execute("""SELECT r.record_id, r.species, r.bin_uri, r.country, r.suspicious_reason
                            FROM records r JOIN sequences s USING(record_id)
                            WHERE r.suspicious=1 AND r.marker_code='COI-5P' AND s.nuc IS NOT NULL""").fetchall()
    n_unmark = 0
    with open(RECHECK, "w", encoding="utf8") as f:
        f.write("record_id\tbin\tspecies\tcountry\tverdict\tchecks now\tdetail now\tmarked because\n")
        for rid, sp, b, co, why in sorted(marked, key=lambda m: m[0]):
            c = cands_all.get(rid, [])
            verdict = "keep" if c else "unmark?"
            n_unmark += not c
            f.write("\t".join([rid, b or "", sp or "", co or "", verdict, ",".join(sorted({x for x, _ in c})),
                               "; ".join(d for _, d in c), (why or "").replace("\t", " ")]) + "\n")
    print(len(marked), "marked records rechecked,", n_unmark, "no longer supported ->", RECHECK)


def main():
    con = sqlite3.connect(br.DB)
    cands, meta = analyse(con, False)
    # A mark counts as supported if either view backs it: judged against the unmarked records only (a bad
    # record cannot excuse itself), or with the marked ones as ordinary records (a name whose records are all
    # marked has no unmarked relatives to compare with).
    cands_strict, meta_all = analyse(con, True, quiet=True, marked_as_evidence=False)
    cands_lenient, _ = analyse(con, True, quiet=True, marked_as_evidence=True)
    cands_all = collections.defaultdict(list)
    for c in (cands_strict, cands_lenient):
        for rid, v in c.items():
            cands_all[rid] += [x for x in v if x not in cands_all[rid]]
    recheck_marked(con, cands_all, meta_all)
    with open(OUT, "w", encoding="utf8") as f:
        f.write("record_id\tbin\tspecies\tcountry\tchecks\tdetail\n")
        for rid in sorted(cands, key=lambda k: (-len({c for c, _ in cands[k]}), k)):
            _, sp, b, co, _ = meta[rid]
            f.write("\t".join([rid, b, sp, co, ",".join(sorted({c for c, _ in cands[rid]})),
                               "; ".join(d for _, d in cands[rid])]) + "\n")
    print(len(cands), "candidates ->", OUT)
    print(collections.Counter(c for v in cands.values() for c in {x for x, _ in v}))


if __name__ == "__main__":
    main()
