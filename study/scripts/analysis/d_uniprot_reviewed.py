"""UniProt reviewed (Swiss-Prot) entries of human, mouse and rat: the lists and the sequences with their modified-residue features.

Query (UniProt REST, endpoint /uniprotkb/stream, format tsv):  (reviewed:true) AND (organism_id:<taxon>)
  rev_<taxon>.tsv  fields accession,reviewed,id,gene_primary   -> columns Entry, Reviewed, Entry Name, Gene Names (primary)
  seq_<taxon>.tsv  fields accession,sequence,ft_mod_res        -> columns Entry, Sequence, Modified residue
Taxa 9606, 10090, 10116. The service serves the current UniProt release; the release number is read from the response header
(X-UniProt-Release) and stored with the SHA-256 of each file. The files used for the stored results are release 2026_03 (see INPUTS.md).

Also writes PRIV/swissprot_reviewed.parquet: the three rev_*.tsv files stacked in the order 9606, 10090, 10116 with a `taxon` column (read by run_amend.py).

Usage: python d_uniprot_reviewed.py [--out DIR] [--no-manifest]
  default: files go to study/analysis/data_private/dl/, the parquet to study/analysis/data_private/, and one record per file is written to
  results/downloads_manifest.json. With --out DIR the files go to DIR and nothing else is written (a test download that leaves the stored files alone).
"""
import argparse
import datetime
import hashlib
import json
import sys

import pandas as pd
import requests

from common import PRIV, RES

URL = "https://rest.uniprot.org/uniprotkb/stream"
TAXA = (9606, 10090, 10116)
KINDS = {"rev": "accession,reviewed,id,gene_primary", "seq": "accession,sequence,ft_mod_res"}


def fetch(kind, taxon, dest):
    params = {"query": f"(reviewed:true) AND (organism_id:{taxon})", "fields": KINDS[kind], "format": "tsv"}
    path = dest / f"{kind}_{taxon}.tsv"
    h = hashlib.sha256()
    with requests.get(URL, params=params, stream=True, timeout=(30, 300)) as r:
        r.raise_for_status()
        head = {k: r.headers.get(k) for k in ("X-UniProt-Release", "X-UniProt-Release-Date", "Date")}
        url = r.url
        with open(path, "wb") as f:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)
                h.update(chunk)
    rows = sum(1 for _ in open(path, "rb")) - 1
    return path, {"source_url": url, "retrieved": datetime.date.today().isoformat(), "uniprot_release": head["X-UniProt-Release"],
                  "uniprot_release_date": head["X-UniProt-Release-Date"], "http_date": head["Date"], "file": path.name, "bytes": path.stat().st_size, "data_rows": rows,
                  "sha256": h.hexdigest()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    ap.add_argument("--no-manifest", action="store_true")
    a = ap.parse_args()
    default = a.out is None
    dest = PRIV / "dl" if default else __import__("pathlib").Path(a.out)
    dest.mkdir(parents=True, exist_ok=True)
    recs = {}
    for kind in KINDS:
        for t in TAXA:
            _, rec = fetch(kind, t, dest)
            recs[f"uniprot_reviewed_{'list' if kind == 'rev' else 'sequences'}_{t}"] = rec
            print(kind, t, rec["uniprot_release"], rec["data_rows"], rec["sha256"], flush=True)
    if default:
        parts = [pd.read_csv(dest / f"rev_{t}.tsv", sep="\t", dtype=str).assign(taxon=t) for t in TAXA]
        pd.concat(parts, ignore_index=True).to_parquet(PRIV / "swissprot_reviewed.parquet")
        if not a.no_manifest:
            man = json.load(open(RES / "downloads_manifest.json"))
            man.update(recs)
            json.dump(man, open(RES / "downloads_manifest.json", "w"), indent=1)
    else:
        json.dump(recs, open(dest / "uniprot_reviewed_download_record.json", "w"), indent=1)


if __name__ == "__main__":
    sys.exit(main())
