"""Build the reference-neighbour file for one specimen from data/bold/bold.sqlite.

Usage: python scripts/build_reference.py [Samlings-nummer|ID] [--k 200] [--out ref]

Files are keyed by the specimen's sequence hash, so specimens sharing a COI sequence share
one file; ref/index.js maps each specimen to its hash (written when building all).

Finds the K reference COI-5P sequences (deduplicated) most similar to the specimen.
Reference sequences have alignment gaps removed and are compared to the specimen at
the best offset within +-MAX_SHIFT bases, since BOLD records start at varying positions.
The output holds, per neighbour: identity, how many records share the sequence, and
their species / BIN / country counts, plus a pairwise distance matrix (tenths of a
percent, specimen = index 0) so the page can run its MDS without any sequences.
"""
import argparse
import collections
import hashlib
import json
import os
import sqlite3
import sys

import numpy as np

DB = "data/bold/bold.sqlite"
MAX_SHIFT = 20
MIN_COMPARED = 400
MAX_LEN = 900
CODE = np.full(256, 4, dtype=np.uint8)
for i, ch in enumerate("ACGT"):
    CODE[ord(ch)] = i
PAD = 5


def load_own():
    s = open("data.js", encoding="utf8").read()
    d = json.loads(s[s.index("["):s.rindex("]") + 1])
    return [r for r in d if r.get("Släkte") == "Coleophora" and r.get("COI")]


def make_slug(all_own):
    """The app's slug rule: collection number when unique, else "<sn>~<ID>", else ID."""
    sn_count = collections.Counter(r["Samlings-nummer"] for r in all_own if r.get("Samlings-nummer"))

    def slug(r):
        sn = r.get("Samlings-nummer")
        return sn if sn and sn_count[sn] == 1 else (sn + "~" + r["ID"] if sn else r["ID"])
    return slug


def encode(seq):
    return CODE[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]


def load_suspicious():
    """The records marked suspicious: (ids, matrix of their sequences). They are kept apart from the
    reference sequences everywhere; pages list the ones near what they show, with the reason."""
    con = sqlite3.connect(DB)
    ids, seqs, meta = [], [], {}
    for rid, nuc, sp, bn in con.execute("""SELECT r.record_id, s.nuc, r.species, r.bin_uri FROM records r JOIN sequences s USING(record_id)
                                   WHERE r.marker_code='COI-5P' AND r.suspicious = 1 AND s.nuc IS NOT NULL"""):
        seq = nuc.upper()
        if MIN_COMPARED <= len(seq.replace("-", "")) and len(seq) <= MAX_LEN:
            ids.append(rid)
            seqs.append(seq)
            meta[rid] = (sp or "", bn or "")
    M = np.full((len(seqs), MAX_LEN), PAD, dtype=np.uint8)
    for i, s in enumerate(seqs):
        M[i, :len(s)] = encode(s)
    return ids, M, meta


def near_suspicious(sus, queries, floor, keep=None, identities=None):
    """Ids of the suspicious sequences at least as alike (best offset) to any of the queries as `floor`;
    keep(species, bin) can narrow them down further. If `identities` is a dict, it gets id -> % identity."""
    ids, M, meta = sus
    if not ids:
        return []
    best = np.zeros(len(ids))
    for q in queries:
        _, ident, _, _ = nearest(q, M, 1)
        best = np.maximum(best, ident)
    found = [ids[i] for i in range(len(ids)) if best[i] >= floor and (keep is None or keep(*meta[ids[i]]))]
    if identities is not None:
        identities.update({i: round(float(best[ids.index(i)]) * 100, 1) for i in found})
    return found


def load_reference():
    """Unique sequences with the records, species, BINs and countries behind each (suspicious records left out)."""
    con = sqlite3.connect(DB)
    groups = collections.OrderedDict()
    q = """SELECT r.record_id, r.species, r.bin_uri, r.country, r.suspicious, s.nuc
           FROM records r JOIN sequences s USING(record_id)
           WHERE r.marker_code='COI-5P' AND s.nuc IS NOT NULL AND r.suspicious = 0"""
    for rid, sp, bn, co, sus, nuc in con.execute(q):
        # Keep BOLD's alignment gaps ("-"): they pad each record to the common frame, so a
        # record starting 25 bases in still lines up. Gaps count as unknown when comparing.
        seq = nuc.upper()
        if not (MIN_COMPARED <= len(seq.replace("-", "")) and len(seq) <= MAX_LEN):
            continue
        g = groups.setdefault(seq, {"ids": [], "species": collections.Counter(),
                                    "bins": collections.Counter(),
                                    "countries": collections.Counter(), "suspicious": 0})
        g["ids"].append(rid)
        g["species"][sp] += 1
        g["bins"][bn] += 1
        g["countries"][co] += 1
        g["suspicious"] += sus
    return list(groups.keys()), list(groups.values())


MAX_INDEL = 3        # longest insertion / deletion (bases) repaired
MIN_INDEL_GAIN = 6   # it must remove at least this many mismatches
MIN_INDEL_TAIL = 15  # ... and be backed by at least this many compared bases after the indel


def find_indel(row, ref):
    """A single insertion/deletion in `row` relative to `ref` (both code arrays in the same frame).

    Past such an indel a column-by-column comparison sees every base as shifted, so one missing
    base near the end of a read looks like dozens of substitutions (a false 5 % distance). Returns
    (p, d) when moving row[p:] by d columns (d > 0: bases missing in row, d < 0: extra bases)
    removes at least MIN_INDEL_GAIN mismatches, else None."""
    L = min(len(row), len(ref))
    row, ref = row[:L], ref[:L]
    ok = (row < 4) & (ref < 4)
    pre = np.concatenate([[0], np.cumsum(ok & (row != ref))])  # pre[i]: mismatches in [0, i)
    base = int(pre[L])
    best = None
    for d in range(-MAX_INDEL, MAX_INDEL + 1):
        if d == 0 or L <= abs(d):
            continue
        a = row[max(0, -d):L - max(0, d)]
        b = ref[max(0, d):L - max(0, -d)]
        okd = (a < 4) & (b < 4)
        mis = okd & (a != b)
        start = max(0, -d)                       # row index of a[0]
        suf = np.zeros(L + 1, dtype=int)         # suf[i]: mismatches of row[i:] against ref[i + d:]
        suf[start:start + len(a)] = np.cumsum(mis[::-1])[::-1]
        cmp_ = np.zeros(L + 1, dtype=int)
        cmp_[start:start + len(a)] = np.cumsum(okd[::-1])[::-1]
        total = pre + suf
        total[cmp_ < MIN_INDEL_TAIL] = base + 1  # too little left after the indel to trust
        i = int(np.argmin(total))
        gain = base - int(total[i])
        if gain >= MIN_INDEL_GAIN and int(suf[i]) <= 0.1 * int(cmp_[i]) and (best is None or gain > best[0]):
            best = (gain, i, d)
    return None if best is None else (best[1], best[2])


def apply_indel(row, p, d, fill):
    """Bring row[p:] back in step: insert d fillers at p (d > 0) or drop -d bases there (d < 0);
    the length stays the same (the end is cut or padded with the filler)."""
    if d > 0:
        return row[:p] + fill * d + row[p:len(row) - d]
    k = -d
    return row[:p] + row[p + k:] + fill * k


def repair_indels(row, ref, fill=None, rounds=2):
    """row repaired against ref (see find_indel); N (code 4) fills the gaps. Works on code arrays or strings."""
    is_str = isinstance(row, str)
    if fill is None:
        fill = "N" if is_str else np.array([4], dtype=np.uint8)
    for _ in range(rounds):
        a = encode(row) if is_str else row
        hit = find_indel(a, ref)
        if hit is None:
            break
        p, d = hit
        row = apply_indel(row, p, d, fill) if is_str else _apply_codes(row, p, d)
    return row


def _apply_codes(row, p, d):
    if d > 0:
        return np.concatenate([row[:p], np.full(d, 4, dtype=row.dtype), row[p:len(row) - d]])
    return np.concatenate([row[:p], row[p - d:], np.full(-d, 4, dtype=row.dtype)])


def nearest(q, R, k):
    """Best-offset identity of query q against every row of R (padded matrix)."""
    n, L = R.shape
    ql = len(q)
    best = np.zeros(n)
    best_cmp = np.zeros(n, dtype=int)
    best_shift = np.zeros(n, dtype=int)
    for s in range(-MAX_SHIFT, MAX_SHIFT + 1):
        # reference position j is compared with query position j + s
        r0, r1 = max(0, -s), min(L, ql - s)
        if r1 <= r0:
            continue
        rs, qs = R[:, r0:r1], q[r0 + s:r1 + s]
        ok = (rs < 4) & (qs < 4)
        compared = ok.sum(axis=1)
        match = ((rs == qs) & ok).sum(axis=1)
        ident = np.where(compared >= MIN_COMPARED, match / np.maximum(compared, 1), 0)
        better = ident > best
        best[better], best_cmp[better], best_shift[better] = ident[better], compared[better], s
    order = np.lexsort((-best_cmp, -best))[:k]
    return order, best, best_cmp, best_shift


def align_to_query(row, shift, ql):
    out = np.full(ql, PAD, dtype=np.uint8)
    r0, r1 = max(0, -shift), min(len(row), ql - shift)
    out[r0 + shift:r1 + shift] = row[r0:r1]
    return out


def distance_matrix(A):
    """Pairwise % difference (x10, int) over positions where both are A/C/G/T."""
    k = A.shape[0]
    D = np.zeros((k, k), dtype=int)
    for i in range(k):
        ok = (A[i] < 4) & (A < 4)
        compared = ok.sum(axis=1)
        diff = ((A[i] != A) & ok).sum(axis=1)
        D[i] = np.where(compared > 0, np.round(1000 * diff / np.maximum(compared, 1)), 1000)
    return D


def build_one(spec_seq, spec, seqs, groups, R, k, sus=None):
    q = encode(spec_seq)
    order, ident, cmp_, shift = nearest(q, R, k)
    A = np.vstack([q] + [repair_indels(align_to_query(R[i, :len(seqs[i])], shift[i], len(q)), q) for i in order])
    D = distance_matrix(A)
    nb = []
    for j, i in enumerate(order):
        g = groups[i]
        nb.append({
            "identity": round(max(float(ident[i]) * 100, 100 - D[0][j + 1] / 10), 2),
            "compared": int(cmp_[i]), "shift": int(shift[i]),
            "n": len(g["ids"]), "ids": g["ids"][:5],
            "species": unnamed(g["species"]), "bins": unnamed(g["bins"]),
            "countries": unnamed(g["countries"]), "suspicious": g["suspicious"],
        })
    out = {"reference_sequences": len(seqs), "k": len(nb), "neighbours": nb, "dist": D.tolist()}
    if sus and nb:
        out["sus_ident"] = {}  # id -> % identity to this specimen
        out["sus"] = near_suspicious(sus, [q], min(float(ident[i]) for i in order), identities=out["sus_ident"])
    return out


def write_js(path, prefix, obj):
    """Plain <script> files (not JSON) so the app also works opened from file://."""
    with open(path, "w", encoding="utf8") as f:
        f.write(prefix + json.dumps(obj, ensure_ascii=False, separators=(",", ":"), sort_keys=True) + ";")


def unnamed(counter):
    """Counter -> dict; missing values (None) get the empty-string key."""
    return {("" if key is None else key): n for key, n in counter.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("specimen", nargs="?", help="Samlings-nummer or ID; default: all")
    ap.add_argument("--k", type=int, default=200)
    ap.add_argument("--out", default="ref")
    a = ap.parse_args()

    own = load_own()
    if a.specimen:
        own = [r for r in own if a.specimen in (r.get("Samlings-nummer"), r.get("ID"))]
        if not own:
            sys.exit("specimen not found (or has no COI): " + a.specimen)

    seqs, groups = load_reference()
    R = np.full((len(seqs), MAX_LEN), PAD, dtype=np.uint8)
    for i, s in enumerate(seqs):
        R[i, :len(s)] = encode(s)
    sus = load_suspicious()

    slug = make_slug(load_own())

    os.makedirs(a.out, exist_ok=True)
    index = {}          # app slug -> sequence hash
    built = {}
    for spec in own:
        qseq = spec["COI"].replace("-", "").upper()
        h = hashlib.sha1(qseq.encode()).hexdigest()[:12]
        index[slug(spec)] = h
        if h in built:
            continue
        built[h] = True
        out = build_one(qseq, spec, seqs, groups, R, a.k, sus)
        out["hash"] = h
        write_js(os.path.join(a.out, h + ".js"), 'window.DNA_REF=window.DNA_REF||{};window.DNA_REF["%s"]=' % h, out)
    if not a.specimen:
        write_js(os.path.join(a.out, "index.js"), "window.DNA_REF_INDEX=", index)
        # The records marked suspicious, with the reason, for the lists of excluded sequences
        con = sqlite3.connect(DB)
        sus = {rid: {"species": sp or "", "country": co or "", "bin": b or "", "reason": why or ""}
               for rid, sp, co, b, why in con.execute(
                   "SELECT record_id, species, country, bin_uri, suspicious_reason FROM records WHERE suspicious = 1")}
        write_js(os.path.join(a.out, "suspicious.js"), "window.DNA_REF_SUS=", sus)
    total = sum(os.path.getsize(os.path.join(a.out, f)) for f in os.listdir(a.out))
    print(len(index), "specimens ->", len(built), "files,", total, "bytes in", a.out)


if __name__ == "__main__":
    main()
