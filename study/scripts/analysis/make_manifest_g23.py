"""Append block G23 (non-primary publications removed from depth; sensitivity only) to grid_manifest.csv: the D0 and D1 rows of G1 and G2 re-run under the
corrected tiers. Added by the plan amendment of 2026-10-04 (sixth batch). Existing rows are not changed. Run once, before any G23 contrast is run."""
import csv

import pandas as pd

from common import REV6

path = REV6 / "grid_manifest.csv"
m = pd.read_csv(path).fillna("")
if (m["block"] == "G23").any():
    raise SystemExit("G23 rows already present")
cols = list(m.columns)
src = m[m["block"].isin(["G1", "G2"]) & m["design"].isin(["D0", "D1"])]
n0 = len(m)
rows = []
for i, (_, r) in enumerate(src.iterrows(), 1):
    x = r.to_dict()
    x.update(block="G23", task="amend", tier_rule="study-cluster rule (884), non-primary PMIDs removed from publication lists",
             variation=f"sensitivity: corrected tiers; source block {r['block']}", status="planned",
             note=f"counterpart of {r['contrast_id']}; added by the plan amendment of 2026-10-04 (sixth batch; scope freeze lifted for G23 only)")
    x["contrast_id"] = f"G23-{n0 + i:04d}"
    rows.append(x)
with open(path, "a", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    for x in rows:
        w.writerow({c: x[c] for c in cols})
print(n0, "->", n0 + len(rows), rows[0]["contrast_id"], rows[-1]["contrast_id"])
