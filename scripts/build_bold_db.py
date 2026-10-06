"""Load BOLD JSON-lines downloads into data/bold/bold.sqlite.

Usage: python scripts/build_bold_db.py [file.jsonl ...]
Defaults to every data/bold/*.jsonl. Re-running upserts records and keeps the
`suspicious` / `suspicious_reason` columns, which are set by hand or by other tools.
"""
import glob
import json
import sqlite3
import sys

DB = "data/bold/bold.sqlite"

SCHEMA = """
CREATE TABLE IF NOT EXISTS records (
    record_id TEXT PRIMARY KEY,          -- processid.marker, unique per sequence
    processid TEXT NOT NULL,
    bin_uri TEXT,
    species TEXT,
    genus TEXT,
    family TEXT,
    identification_rank TEXT,
    country TEXT,
    province TEXT,
    lat REAL,
    lon REAL,
    collection_date_start TEXT,
    marker_code TEXT,
    nuc_basecount INTEGER,
    suspicious INTEGER NOT NULL DEFAULT 0,
    suspicious_reason TEXT
);
CREATE TABLE IF NOT EXISTS sequences (
    record_id TEXT PRIMARY KEY REFERENCES records(record_id),
    nuc TEXT
);
CREATE TABLE IF NOT EXISTS raw (
    record_id TEXT PRIMARY KEY REFERENCES records(record_id),
    json TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_records_bin ON records(bin_uri);
CREATE INDEX IF NOT EXISTS idx_records_species ON records(species);
CREATE INDEX IF NOT EXISTS idx_records_country ON records(country);
CREATE INDEX IF NOT EXISTS idx_records_marker ON records(marker_code);
CREATE INDEX IF NOT EXISTS idx_records_processid ON records(processid);

-- Which BINs does a species occur in, and vice versa (named records only).
CREATE VIEW IF NOT EXISTS species_bins AS
    SELECT species, bin_uri, COUNT(*) AS n FROM records
    WHERE species IS NOT NULL AND bin_uri IS NOT NULL
    GROUP BY species, bin_uri;

-- Every record in any BIN that contains the species, including unnamed and
-- differently named records. Join on species_name to get the comparison set.
CREATE VIEW IF NOT EXISTS species_comparison_set AS
    SELECT sb.species AS species_name, r.*
    FROM (SELECT DISTINCT species, bin_uri FROM species_bins) sb
    JOIN records r ON r.bin_uri = sb.bin_uri;
"""

UPSERT = """
INSERT INTO records (record_id, processid, bin_uri, species, genus, family,
    identification_rank, country, province, lat, lon, collection_date_start,
    marker_code, nuc_basecount)
VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)
ON CONFLICT(record_id) DO UPDATE SET
    processid=excluded.processid, bin_uri=excluded.bin_uri, species=excluded.species,
    genus=excluded.genus, family=excluded.family,
    identification_rank=excluded.identification_rank, country=excluded.country,
    province=excluded.province, lat=excluded.lat, lon=excluded.lon,
    collection_date_start=excluded.collection_date_start,
    marker_code=excluded.marker_code, nuc_basecount=excluded.nuc_basecount
"""


def blank(v):
    return v if v not in ("", None) else None


def load(con, path):
    n = 0
    with open(path, encoding="utf8") as f:
        for line in f:
            r = json.loads(line)
            rid = r["record_id"]
            coord = r.get("coord") or [None, None]
            con.execute(UPSERT, (
                rid, r["processid"], blank(r.get("bin_uri")), blank(r.get("species")),
                blank(r.get("genus")), blank(r.get("family")),
                blank(r.get("identification_rank")), blank(r.get("country/ocean")),
                blank(r.get("province/state")), coord[0], coord[1],
                blank(r.get("collection_date_start")), blank(r.get("marker_code")),
                r.get("nuc_basecount"),
            ))
            con.execute("INSERT OR REPLACE INTO sequences VALUES (?,?)",
                        (rid, blank(r.get("nuc"))))
            con.execute("INSERT OR REPLACE INTO raw VALUES (?,?)", (rid, line.strip()))
            n += 1
    return n


def main():
    files = sys.argv[1:] or sorted(glob.glob("data/bold/*.jsonl"))
    con = sqlite3.connect(DB)
    con.executescript(SCHEMA)
    for path in files:
        print(path, load(con, path))
    con.commit()
    con.execute("ANALYZE")
    con.close()


if __name__ == "__main__":
    main()
