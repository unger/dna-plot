# Coleophora DNA

A small static web app that lists all *Coleophora* records from the DNA dataset,
with a free-text filter and a detail page showing the full data including the COI
sequence and source link.

**Live site:** https://unger.github.io/dna-plot/

## Files

- `index.html` — the whole app (HTML + CSS + JS, no dependencies).
- `data.js` — the dataset exposed as `window.DNA_DATA`. Loaded via a `<script>`
  tag so the app works both when opened directly (`file://`) and when hosted.
- `doc/Dna svar.json` — the source data (`data.js` is generated from this).

## Run locally

Just open `index.html` in a browser — no server needed.

## Regenerate `data.js` after the data changes

```bash
python -c "import json; d=json.load(open('doc/Dna svar.json',encoding='utf-8')); \
open('data.js','w',encoding='utf-8').write('window.DNA_DATA='+json.dumps(d,ensure_ascii=False,separators=(',',':'))+';')"
```

## Host on GitHub Pages

1. Create a repository on GitHub and push these files:
   ```bash
   git init
   git add .
   git commit -m "Coleophora DNA app"
   git branch -M main
   git remote add origin https://github.com/<user>/<repo>.git
   git push -u origin main
   ```
2. On GitHub: **Settings → Pages → Build and deployment**, set
   **Source = Deploy from a branch**, **Branch = `main`**, folder **`/ (root)`**, Save.
3. The site appears at `https://<user>.github.io/<repo>/` within a minute or two.

The `.nojekyll` file tells GitHub Pages to serve the files as-is.

## BIN reference page

`#/b` lets you pick one or several BOLD BINs (at most 8) and plots only their COI sequences,
one color per BIN, with your own findings as black diamonds and a heading per BIN listing the
species in it. Your findings are placed in the BIN of their nearest BOLD sequence (cached in
`ref/summary.js`). Data is built from `data/bold/bold.sqlite`:

```bash
python scripts/build_bin_reference.py   # after build_ref_summary.py and build_species_reference.py
```

## Data source

Reference sequences, species names and BIN assignments come from BOLD Systems
(https://boldsystems.org), retrieved 2026-10-06 (genus *Coleophora*, marker COI-5P).
Please cite, as BOLD asks (https://boldsystems.org/about/citation/):

- Ratnasingham S, Hebert PDN (2013). A DNA-Based Registry for All Animal Species: The Barcode Index
  Number (BIN) System. PLoS ONE 8(8): e66213. doi:10.1371/journal.pone.0066213
- Ratnasingham S, Wei C, Chan D, et al., Hebert PDN (2024). BOLD v4: A Centralized Bioinformatics Platform
  for DNA-Based Biodiversity Data. In: DNA Barcoding: Methods and Protocols, Springer, pp. 403-441.

Individual records have no DOI; cite them as "ProcessID. Specimen Depository. (Accessed on YYYY-MM-DD) via
boldsystems.org". The BIN of the own findings is an estimate from the nearest BOLD sequence within 2.5 %
difference, not BOLD's own assignment.
