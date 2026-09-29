# Coleophora DNA

A small static web app that lists all *Coleophora* records from the DNA dataset,
with a free-text filter and a detail page showing the full data including the COI
sequence and source link.

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
