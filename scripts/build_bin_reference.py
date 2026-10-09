"""Reference data per BOLD BIN: ref/bin_<id>.js (the BIN's sequences) and ref/bins.js (overview).

Usage: python scripts/build_bin_reference.py   (after build_ref_summary.py and build_species_reference.py)

Covers the BINs the own individuals fall in (BIN of the nearest BOLD sequence, cached in
ref/summary.js), the BINs of the own species names (ref/species.js) and the BINs of the species
listed in scripts/bin_species.txt (one name per line), so a species you have no findings of can be
looked up too. Per BIN file:
  * seqs - the BIN's unique COI sequences in BOLD's alignment frame (gaps "-"), each with the
    records, species, countries and ids behind it
  * own  - the own individuals assigned to the BIN, put in the same frame, identical ones merged
The page picks one or more BINs and computes the distances between all these sequences itself.
"""
import collections
import json
import sqlite3

import numpy as np

import build_reference as br

TOP = 4  # species listed per BIN in the overview
NEAR = 25  # nearest other BINs listed per BIN
FRAME = 700  # columns of the alignment frame the consensus sequences cover
MAX_FRAME_SHIFT = 450  # BOLD does not always pad a record that starts late; look this far for its place
MIN_FRAME_IDENT = 0.9  # a record this alike the anchor as it stands is in the frame
MIN_FRAME_COMPARED = 120


def own_in_frame(qseq, R):
    """The own sequence (no gaps) placed in BOLD's alignment frame, via the best-matching reference."""
    q = br.encode(qseq)
    order, _, _, shift = br.nearest(q, R, 1)
    s = int(shift[order[0]])
    # reference position j is compared with query position j + s, so q[p] sits at frame column p - s
    start = max(0, -s)
    return "-" * start + qseq[max(0, s):]


def frame_identity(q, anchor, k):
    """Identity of q placed at frame column p + k against the anchor, and the number of compared bases."""
    lo, hi = max(0, -k), min(len(q), len(anchor) - k)
    if hi <= lo:
        return 0.0, 0
    a, b = q[lo:hi], anchor[lo + k:hi + k]
    ok = (a < 4) & (b < 4)
    n = int(ok.sum())
    return (float(((a == b) & ok).sum()) / n if n else 0.0), n


def in_bin_frame(seq, anchor):
    """Put a record where it lines up with the BIN's anchor sequence. Most records already sit in BOLD's
    frame, but some start late without padding, and some have BOLD's gaps in the wrong places; those are
    stripped of gaps and placed at the offset where they match the anchor best."""
    q = br.encode(seq)  # keeps "-" (code 4): leading gaps already place the record in the frame
    ident0, _ = frame_identity(q, anchor, 0)
    if ident0 >= MIN_FRAME_IDENT:
        return seq
    stripped = seq.replace("-", "")
    qs = br.encode(stripped)
    lead = len(seq) - len(seq.lstrip("-"))
    best, best_k = -1.0, 0
    for k in [lead] + [k for k in range(-MAX_FRAME_SHIFT, MAX_FRAME_SHIFT + 1) if k != lead]:
        ident, n = frame_identity(qs, anchor, k)
        if n >= MIN_FRAME_COMPARED and ident > best:
            best, best_k = ident, k
        if k == lead and best >= MIN_FRAME_IDENT:
            break
    if best < ident0 + 0.05:
        return seq
    return "-" * best_k + stripped if best_k > 0 else stripped[-best_k:]


def consensus(by_seq):
    """The commonest base per frame column over the BIN's records (4 = nothing there)."""
    counts = np.zeros((FRAME, 4))
    for seq, g in by_seq.items():
        codes = br.encode(seq)[:FRAME]
        for c in range(4):
            counts[:len(codes), c] += (codes == c) * len(g["ids"])
    return np.where(counts.sum(axis=1) > 0, counts.argmax(axis=1), 4).astype(np.uint8)


def nearest_bins(cons):
    """For each BIN the NEAR closest others as [bin, tenths of a percent] between consensus sequences."""
    names = list(cons)
    M = np.vstack([cons[b] for b in names])
    near = {}
    for i, b in enumerate(names):
        ok = (M[i] < 4) & (M < 4)
        compared = ok.sum(axis=1)
        diff = ((M[i] != M) & ok).sum(axis=1)
        d = np.where(compared >= br.MIN_COMPARED // 2, np.round(1000 * diff / np.maximum(compared, 1)), 10000)
        order = [j for j in np.argsort(d, kind="stable") if j != i and d[j] < 10000][:NEAR]
        near[b] = [[names[j], int(d[j])] for j in order]
    return near


def main():
    own = br.load_own()
    slug = br.make_slug(own)
    summary = json.loads(open("ref/summary.js", encoding="utf8").read().split("=", 1)[1].rstrip(";"))
    species = json.loads(open("ref/species.js", encoding="utf8").read().split("=", 1)[1].rstrip(";"))

    own_by_bin = collections.defaultdict(lambda: collections.OrderedDict())  # bin -> sequence -> [slugs]
    for rec in own:
        b = summary.get(slug(rec), {}).get("bin")
        if b:
            own_by_bin[b].setdefault(rec["COI"].replace("-", "").upper(), []).append(slug(rec))
    wanted = set(own_by_bin) | {b["bin"] for sp in species.values() for b in sp["bins"]}
    con = sqlite3.connect(br.DB)
    for line in open("scripts/bin_species.txt", encoding="utf8"):
        name = line.strip()
        if name and not name.startswith("#"):
            wanted |= {r[0] for r in con.execute(
                "SELECT DISTINCT bin_uri FROM records WHERE species = ? AND bin_uri IS NOT NULL", (name,))}

    seqs, groups = br.load_reference()
    R = np.full((len(seqs), br.MAX_LEN), br.PAD, dtype=np.uint8)
    for i, s in enumerate(seqs):
        R[i, :len(s)] = br.encode(s)

    overview = {}
    cons = {}
    for b in sorted(wanted):
        rows = con.execute(
            """SELECT r.record_id, r.species, r.country, r.suspicious, s.nuc, r.suspicious_reason FROM records r
               JOIN sequences s USING(record_id)
               WHERE r.bin_uri = ? AND r.marker_code = 'COI-5P' AND s.nuc IS NOT NULL""", (b,)).fetchall()
        flagged = [{"id": rid, "species": sp or "", "country": co or "", "reason": why or ""}
                   for rid, sp, co, sus, nuc, why in rows if sus]  # listed on the BIN page with the reason
        rows = [r for r in rows if not r[3]]  # suspicious records stay out of everything below
        anchor_seq = collections.Counter(r[4].upper() for r in rows if "-" not in r[4] and len(r[4]) >= 640).most_common(1) or \
            collections.Counter(r[4].upper() for r in rows).most_common(1)
        anchor = br.encode(anchor_seq[0][0]) if anchor_seq else None
        by_seq = collections.OrderedDict()
        for rid, sp, co, sus, nuc, why in rows:
            nuc = nuc.upper()
            if anchor is not None:
                nuc = br.repair_indels(in_bin_frame(nuc, anchor), anchor)
            g = by_seq.setdefault(nuc, {"ids": [], "species": collections.Counter(),
                                        "countries": collections.Counter(), "suspicious": 0})
            g["ids"].append(rid)
            g["species"][sp] += 1
            g["countries"][co] += 1
        out_seqs = [{"s": s, "n": len(g["ids"]), "ids": g["ids"][:5], "species": br.unnamed(g["species"]),
                     "countries": br.unnamed(g["countries"]), "suspicious": g["suspicious"]}
                    for s, g in sorted(by_seq.items(), key=lambda kv: -len(kv[1]["ids"]))]
        out_own = [{"s": br.repair_indels(in_bin_frame(own_in_frame(q, R), anchor), anchor) if anchor is not None else own_in_frame(q, R), "slugs": sl} for q, sl in own_by_bin.get(b, {}).items()]
        key = b.split(":")[-1]
        br.write_js("ref/bin_%s.js" % key, 'window.DNA_REF_BIN=window.DNA_REF_BIN||{};window.DNA_REF_BIN["%s"]=' % b,
                    {"bin": b, "seqs": out_seqs, "own": out_own, "sus": flagged})

        cons[b] = consensus(by_seq)
        sp_total = collections.Counter()
        for g in by_seq.values():
            sp_total.update(g["species"])
        overview[b] = {"file": key, "n": len(rows), "variants": len(out_seqs),
                       "species": br.unnamed(collections.Counter(dict(sp_total.most_common(TOP)))),
                       "n_species": len(sp_total),
                       "names": [k for k, _ in sp_total.most_common() if k], "own": sum(len(v) for v in own_by_bin.get(b, {}).values())}
        print(b, len(rows), "records,", len(out_seqs), "sequences,", overview[b]["own"], "own")
    for b, near in nearest_bins(cons).items():
        overview[b]["near"] = near
    br.write_js("ref/bins.js", "window.DNA_REF_BINS=", overview)


if __name__ == "__main__":
    main()
