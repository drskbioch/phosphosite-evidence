"""G23 step 1: candidate PMIDs (study-cluster candidates, linked to >= 884 records) plus PMID 25514926, and their PubMed metadata (esummary).
Output results/g23_pubmed_metadata.json. No depth or tier is recomputed here."""
import json
import time
import urllib.parse
import urllib.request
from datetime import date

import pandas as pd

from common import REL, RES, cluster_map, split_pmids

m = pd.read_parquet(REL / "mammalian_sites.parquet", columns=["pmid_list_str"])
L = split_pmids(m["pmid_list_str"])
cid, info, freq = cluster_map(L, 884)
cand = sorted(cid)
extra = "25514926"
cand_all = sorted(set(cand) | {extra})
members = {}
for p, c in cid.items():
    members.setdefault(c, []).append(p)
out = {"retrieved": str(date.today()), "n_candidates_study_cluster_rule": len(cand), "extra_pmid_in_candidates": extra in cid, "candidates": cand_all,
       "records_linked": {p: int(freq[p]) for p in cand_all}, "cluster_of_25514926": sorted(members[cid[extra]]) if extra in cid else None, "records": {}, "missing": []}
URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
for i in range(0, len(cand_all), 100):
    batch = cand_all[i:i + 100]
    data = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(batch), "retmode": "json"}).encode()
    for attempt in range(4):
        try:
            r = json.load(urllib.request.urlopen(urllib.request.Request(URL, data=data), timeout=60))
            break
        except Exception as e:  # noqa: BLE001
            print("retry", i, e)
            time.sleep(3)
    res = r["result"]
    for p in batch:
        d = res.get(p)
        if not d or "error" in d:
            out["missing"].append(p)
            continue
        out["records"][p] = {"title": d.get("title"), "source": d.get("source"), "pubdate": d.get("pubdate"), "pubtype": d.get("pubtype"), "elocationid": d.get("elocationid")}
    time.sleep(0.5)
json.dump(out, open(RES / "g23_pubmed_metadata.json", "w"), indent=1)
print("candidates", len(cand), "plus extra in set:", out["extra_pmid_in_candidates"], "| fetched", len(out["records"]), "| missing", out["missing"])
print("cluster of 25514926:", out["cluster_of_25514926"])
