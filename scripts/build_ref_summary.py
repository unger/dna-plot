"""Per-individual check of the species name against the BOLD reference: ref/summary.js.

Usage: python scripts/build_ref_summary.py   (after build_reference.py has written ref/)

For each individual, compares its scientific name with the species names of its nearest
reference sequences. Status:
  ok     the individual's name is among the names of the nearest named sequences
  stav   the same name spelled differently in BOLD (doubled letters, e.g. galbulipennella)
  annat  another name is clearly nearer (or the individual's own name is not nearby at all)
  langt  no named sequence within FAR % difference, so nothing to judge by
  okand  the individual has no species epithet to compare
Also gives the BIN of the nearest sequence that has one.
"""
import json
import re

from build_reference import load_own, make_slug, write_js

TOLERANCE = 0.3   # percentage points: sequences this close to the best count as equally near
FAR = 2.5         # nearest named sequence beyond this % difference = nothing to judge by
BIN_MAX = 2.5     # only assign a probable BIN when the nearest sequence with one is this close


def load_js(path):
    s = open(path, encoding="utf8").read()
    return json.loads(s[s.index("={") + 1:s.rindex(";")])


def norm(name):
    return re.sub(r"\s+", " ", (name or "").replace("?", "")).strip()


def fold(name):
    """Spelling-tolerant key: lower case, repeated letters collapsed."""
    return re.sub(r"(.)\1+", r"\1", norm(name).lower())


def diff(nb):
    return round(100 - nb["identity"], 2)


def main():
    index = load_js("ref/index.js")
    slug = make_slug(load_own())
    summary, cache = {}, {}
    for rec in load_own():
        sl = slug(rec)
        h = index.get(sl)
        if not h:
            continue
        if h not in cache:
            # sequences only suspicious records have are not evidence for or against a name
            cache[h] = [nb for nb in load_js("ref/%s.js" % h)["neighbours"] if nb["suspicious"] < nb["n"]]
        nbs = cache[h]
        own = norm(rec.get("Vetenskapligt namn"))
        named = [nb for nb in nbs if any(k for k in nb["species"])]
        entry = {"h": h}
        if named:
            d0 = diff(named[0])
            near = [nb for nb in named if diff(nb) <= d0 + TOLERANCE]
            best = {}
            for nb in near:
                for k, n in nb["species"].items():
                    if k:
                        best[k] = best.get(k, 0) + n
            own_d = next((diff(nb) for nb in named if own in nb["species"]), None)
            folded = {fold(k): k for k in best}
            entry.update(d0=d0, best=best, own_d=own_d)
        if len(own.split()) < 2:
            entry["st"] = "okand"
        elif not named or entry["d0"] > FAR:
            entry["st"] = "langt"
        elif own in entry["best"]:
            entry["st"] = "ok"
        elif fold(own) in folded:
            entry["st"] = "stav"
        else:
            entry["st"] = "annat"
        with_bin = next((nb for nb in nbs if any(k for k in nb["bins"])), None)
        if with_bin and diff(with_bin) <= BIN_MAX:
            top = max((k for k in with_bin["bins"] if k), key=lambda k: with_bin["bins"][k])
            entry["bin"] = top
            entry["bin_d"] = diff(with_bin)
        summary[sl] = entry
    write_js("ref/summary.js", "window.DNA_REF_SUMMARY=", summary)
    counts = {}
    for e in summary.values():
        counts[e["st"]] = counts.get(e["st"], 0) + 1
    print(len(summary), "individuals;", counts)


if __name__ == "__main__":
    main()
