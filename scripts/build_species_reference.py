"""Reference data per species: ref/sp_<hash>.js (plot data) and ref/species.js (BIN overview).

Usage: python scripts/build_species_reference.py   (after build_reference.py and build_ref_summary.py)

Per species name in the app:
  * ref/species.js  - how the name is used in BOLD: records, BINs, other names sharing those BINs
  * ref/sp_<hash>.js - the species' own COI variants plus up to K_TOTAL reference sequences
    (those in the species' BINs first, then the nearest others) and their pairwise distances,
    in the same format as the per-individual files, with `own` listing the species' variants.
"""
import hashlib
import json
import re
import sqlite3
import collections

import numpy as np

import build_reference as br

K_TOTAL = 150     # reference sequences per species file
K_BIN = 100       # at most this many of them chosen because they sit in the species' BINs
TOP = 6           # other names / countries listed per BIN


def norm(name):
    return re.sub(r"\s+", " ", (name or "").replace("?", "")).strip()


def top(counter, n):
    return dict(counter.most_common(n))


def bin_overview(con, name):
    """How BOLD uses the species name: records and BINs, and who else is in those BINs."""
    rows = con.execute("SELECT bin_uri, COUNT(*) FROM records WHERE species = ? GROUP BY bin_uri", (name,)).fetchall()
    bins = []
    for bin_uri, n in sorted(rows, key=lambda r: -r[1]):
        if bin_uri is None:
            continue
        others = collections.Counter()
        countries = collections.Counter()
        total = 0
        for sp, co, c in con.execute(
                "SELECT species, country, COUNT(*) FROM records WHERE bin_uri = ? GROUP BY species, country", (bin_uri,)):
            total += c
            countries[co or ""] += c
            if sp != name:
                others[sp or ""] += c
        bins.append({"bin": bin_uri, "n": n, "total": total,
                     "others": top(others, TOP), "n_others": len(others), "countries": top(countries, TOP)})
    return {"named": sum(r[1] for r in rows), "no_bin": sum(n for b, n in rows if b is None), "bins": bins}


def main():
    own = br.load_own()
    slug = br.make_slug(own)
    summary = json.loads(open("ref/summary.js", encoding="utf8").read().split("=", 1)[1].rstrip(";"))
    con = sqlite3.connect(br.DB)
    seqs, groups = br.load_reference()
    R = np.full((len(seqs), br.MAX_LEN), br.PAD, dtype=np.uint8)
    for i, s in enumerate(seqs):
        R[i, :len(s)] = br.encode(s)

    by_name = collections.OrderedDict()
    for rec in own:
        by_name.setdefault(rec.get("Vetenskapligt namn") or "Okänd art", []).append(rec)

    overview = {}
    for name, members in by_name.items():
        bold_name = norm(name)
        if len(bold_name.split()) < 2:
            continue                       # genus only: nothing to compare by name
        h = hashlib.sha1(name.encode("utf8")).hexdigest()[:10]

        # The species' own variants, identical sequences merged
        variants = collections.OrderedDict()
        for r in members:
            variants.setdefault(r["COI"].replace("-", "").upper(), []).append(slug(r))
        vseqs = list(variants.keys())
        q = [br.encode(v) for v in vseqs]

        # BINs the species is in: from the name in BOLD and from the individuals' nearest sequences
        ov = bin_overview(con, bold_name)
        species_bins = {b["bin"] for b in ov["bins"]}
        for sl in sum(variants.values(), []):
            b = summary.get(sl, {}).get("bin")
            if b:
                species_bins.add(b)
        ov["file"] = h
        ov["own_bins"] = collections.Counter(summary.get(sl, {}).get("bin") for sl in sum(variants.values(), []))
        ov["own_bins"] = {(k or ""): v for k, v in ov["own_bins"].items()}
        overview[name] = ov

        # Best identity of each reference sequence to any of the variants
        best = np.zeros(len(seqs))
        shift0 = None
        for qi, qv in enumerate(q):
            _, ident, _, shift = br.nearest(qv, R, 1)
            best = np.maximum(best, ident)
            if qi == 0:
                shift0 = shift                  # frame: the first variant
        in_bin = np.array([any(b in g["bins"] for b in species_bins) for g in groups])
        order_all = np.argsort(-best, kind="stable")
        chosen = [i for i in order_all if in_bin[i]][:K_BIN]
        seen = set(chosen)
        chosen += [i for i in order_all if i not in seen][:max(0, K_TOTAL - len(chosen))]
        chosen.sort(key=lambda i: -best[i])

        # Align everything to the first variant's frame, then pairwise distances
        ql = len(vseqs[0])
        Rown = np.full((len(q), br.MAX_LEN), br.PAD, dtype=np.uint8)
        for i, qv in enumerate(q):
            Rown[i, :len(qv)] = qv
        _, _, _, own_shift = br.nearest(q[0], Rown, len(q))
        rows = [br.align_to_query(Rown[i, :len(q[i])], own_shift[i], ql) for i in range(len(q))]
        rows += [br.align_to_query(R[i, :len(seqs[i])], shift0[i], ql) for i in chosen]
        D = br.distance_matrix(np.vstack(rows))

        nb = []
        for i in chosen:
            g = groups[i]
            nb.append({
                "identity": round(float(best[i]) * 100, 2), "n": len(g["ids"]), "ids": g["ids"][:5],
                "species": br.unnamed(g["species"]), "bins": br.unnamed(g["bins"]),
                "countries": br.unnamed(g["countries"]), "suspicious": g["suspicious"],
                "inbin": bool(in_bin[i]),
            })
        out = {"own": [{"slugs": v, "n": len(v)} for v in variants.values()],
               "reference_sequences": len(seqs), "k": len(nb), "neighbours": nb, "dist": D.tolist()}
        br.write_js("ref/sp_%s.js" % h, 'window.DNA_REF_SP=window.DNA_REF_SP||{};window.DNA_REF_SP["%s"]=' % h, out)
        print(name, len(variants), "variants,", len(nb), "refs,", int(in_bin[chosen].sum()), "in BINs")
    br.write_js("ref/species.js", "window.DNA_REF_SPECIES=", overview)


if __name__ == "__main__":
    main()
