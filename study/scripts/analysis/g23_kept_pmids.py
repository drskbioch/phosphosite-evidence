"""G23 addendum (descriptive only; the G23 rule and every tier are unchanged): records whose publication list contains each of three PMIDs that were kept as primary
although their titles are not phosphorylation evidence: 30243273, 37031204, 14702039. Output results/g23_kept_nonevidence_pmids.json."""
import json

import numpy as np
import pandas as pd

from common import REL, RES

meta = json.load(open(RES / "g23_pubmed_metadata.json"))
m = pd.read_parquet(REL / "mammalian_sites.parquet", columns=["pmid_list_str", "tier", "species"])
s = m["pmid_list_str"].fillna("")
out = {}
for p in ("30243273", "37031204", "14702039"):
    has = s.str.contains(rf"(?:^|;){p}(?:;|$)", regex=True).to_numpy()
    out[p] = {"title": meta["records"][p]["title"], "records": int(has.sum()), "share_of_records": float(has.mean()), "by_tier": pd.Series(m["tier"].to_numpy()[has]).value_counts().to_dict(),
              "by_species": pd.Series(m["species"].to_numpy()[has]).value_counts().to_dict(), "only_pmid_in_list": int((s.to_numpy()[has] == p).sum())}
any_ = np.zeros(len(m), bool)
for p in out:
    any_ |= s.str.contains(rf"(?:^|;){p}(?:;|$)", regex=True).to_numpy()
out["any_of_the_three"] = {"records": int(any_.sum()), "share_of_records": float(any_.mean())}
json.dump(out, open(RES / "g23_kept_nonevidence_pmids.json", "w"), indent=1, default=int)
print(json.dumps(out, indent=1)[:1800])
