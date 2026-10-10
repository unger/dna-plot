"""Suggest suspicious BOLD records for review: data/bold/suspicious_candidates.tsv.

Usage: python scripts/find_suspicious.py   (reads data/bold/bold.sqlite, changes nothing)

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
NAME_MARGIN = 1.5   # % the nearest other-species record must be closer than the nearest same-species one
MAJORITY_SHARE = 0.6    # a BIN's majority name is at least this share of its named records
MINORITY_SHARE = 0.2    # bin_minority_name: the record's name is at most this share
SWAP_SHARE = 0.5        # name_swap: looser, for small BINs (e.g. 1 of 3)
ESTABLISHED = 3         # name_outlier: a name with this many records in a BIN, at least half of its named ones, is not a slip
MIN_COMPARED = 400
OUT = "data/bold/suspicious_candidates.tsv"
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


def main():
    con = sqlite3.connect(br.DB)
    rows = con.execute(
        """SELECT r.record_id, r.species, r.bin_uri, r.country, s.nuc FROM records r
           JOIN sequences s USING(record_id)
           WHERE r.marker_code='COI-5P' AND s.nuc IS NOT NULL AND r.suspicious=0""").fetchall()
    rows = [(rid, sp or "", b or "", co or "", nuc.upper()) for rid, sp, b, co, nuc in rows if len(nuc) <= br.MAX_LEN]
    fixes = name_fixes(rows)
    for old, (new, why) in sorted(fixes.items()):
        print("%s -> %s (%s, %d record(s))" % (old, new, why, sum(1 for r in rows if r[1] == old)))
    rows = [(rid, fixes.get(sp, (sp,))[0], b, co, nuc) for rid, sp, b, co, nuc in rows]
    cands = collections.defaultdict(list)  # record_id -> [(check, detail)]
    meta = {r[0]: r for r in rows}

    # --- stop codons
    for rid, n in stop_candidates(rows).items():
        cands[rid].append(("stop_codon", "%d stop codon(s) in the best reading frame" % n))

    # --- put every record in its BIN's frame (BOLD does not pad all records), so columns can be compared
    def anchor_of(seqs):
        c = collections.Counter(x for x in seqs if "-" not in x and len(x) >= 640) or collections.Counter(seqs)
        return br.encode(c.most_common(1)[0][0])
    by_bin_rows = collections.defaultdict(list)
    for r in rows:
        by_bin_rows[r[2]].append(r)
    glob = anchor_of([r[4] for r in rows])
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

    # --- BIN outliers: nearest other record in the same BIN
    by_bin = collections.defaultdict(set)
    for r in rows:
        if r[2]:
            by_bin[r[2]].add(idx[r[4]])
    for b, members in by_bin.items():
        members = sorted(members)
        if len(members) < 4:
            continue
        for i in members:
            d = pdist(R[i], R[[j for j in members if j != i]])
            others = [j for j in members if j != i]
            order = np.argsort(np.where(np.isnan(d), 1e9, d))[:5]
            best = np.nanmin(d) if not np.all(np.isnan(d)) else np.nan
            if np.isnan(best) or best <= BIN_OUTLIER:
                continue
            best = min(repaired_dist(R[i], R[others[k]]) if not np.isnan(d[k]) else 1e9 for k in order)
            if best > BIN_OUTLIER:
                for r in recs_of[i]:
                    cands[r[0]].append(("bin_outlier", "%.1f %% from the nearest of %d others in %s" % (best, len(members) - 1, b)))

    # --- name outliers: nearest same-species vs nearest other-species record
    n_name_bin = collections.Counter((r[2], r[1]) for r in rows if r[2] and proper_name(r[1]))
    n_bin = collections.Counter(r[2] for r in rows if r[2] and proper_name(r[1]))

    def established(r):
        """Several records with the same name in the BIN, at least half of its named records: more than a slip,
        so the other-named records there are the odd ones (or a species COI barely separates)."""
        return n_name_bin[r[2], r[1]] >= ESTABLISHED and 2 * n_name_bin[r[2], r[1]] >= n_bin[r[2]]

    mask = collections.defaultdict(lambda: np.zeros(len(useqs), dtype=bool))
    for j, nm in enumerate(names_of):
        for n in nm:
            mask[n][j] = True
    named = np.array([bool(nm) for nm in names_of])
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
                near = ", ".join(sorted(names_of[jo[int(np.nanargmin([repaired_dist(R[i], R[j]) for j in jo]))]]))
                own = "%.1f %% from the nearest of its own name" % rs
                for r in recs_of[i]:
                    if established(r):
                        continue
                    cands[r[0]].append(("name_outlier", "%s; %.1f %% from %s" % (own, ro, near)))
    # --- a name in a BIN of another species, and exchanged names. Judged per record rather than per unique
    # sequence, so a sequence shared by several species does not hide a record's own name.
    seqs_of = collections.defaultdict(lambda: collections.defaultdict(set))  # name -> bin -> unique sequences
    named_in = collections.defaultdict(list)  # bin -> records with a proper name
    for r in rows:
        if proper_name(r[1]):
            seqs_of[r[1]][r[2]].add(idx[r[4]])
            if r[2]:
                named_in[r[2]].append(r)
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
    for b, mem in named_in.items():
        cnt = collections.Counter(r[1] for r in mem)
        top = cnt.most_common(2)
        M, mc = top[0]
        if len(mem) < 3 or (len(top) > 1 and top[1][1] == mc) or mc / len(mem) < MAJORITY_SHARE:
            continue
        for r in mem:
            X = r[1]
            if X == M or cnt[X] >= mc or cnt[X] / len(mem) > SWAP_SHARE:
                continue
            i = idx[r[4]]
            elsewhere = set().union(*[v for k, v in seqs_of[X].items() if k != b])
            dM = nearest(i, seqs_of[M][b], ("in", M, b))
            dX = nearest(i, elsewhere, ("out", X, b))
            if dM <= BIN_OUTLIER and np.isfinite(dX) and dX - dM >= NAME_MARGIN:
                dissent.append((r, X, M, b, mc, len(mem), cnt[X], dM, dX))
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
