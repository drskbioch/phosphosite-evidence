"""Append blocks G21 (within-gene GOLD vs BRONZE) and G22 (post-hoc S/T-only sequence-motif contrast) to grid_manifest.csv.
Added by the plan amendment of 2026-10-04 after the project lead reviewed RESULTS_REPORT.md. Existing rows are not changed. Run once, before any contrast
of these blocks is run. Idempotent: refuses if G21 rows already exist."""
import csv

import pandas as pd

from common import REV6

path = REV6 / "grid_manifest.csv"
m = pd.read_csv(path).fillna("")
if (m["block"] == "G21").any():
    raise SystemExit("G21 rows already present")
cols = list(m.columns)
SIG = {"S2": "conservation", "S3": "disorder", "S4": "kinase motif (sequence)", "S5": "kinase motif (lookup)"}
NOTE = "added by the plan amendment of 2026-10-04 (after the results report was reviewed)"
rows = []
MODELS = [("M1", "M1: all records of the population"), ("M1s", "M1s: records with a protein sequence, no detectability control"),
          ("M2", "M2: records with a protein sequence, with the detectability control")]
for pop in ("P-all (S/T/Y)", "P-human (S/T/Y)"):
    for s in SIG:
        for key, lab in MODELS:
            rows.append(dict(block="G21", task="amend", signal=s, signal_label=SIG[s], population=pop, design="D6 (within-gene)", index="not used",
                             tier_rule="study-cluster rule (884)", variation=lab, decision_bearing="no" if s == "S5" else "yes", status="planned",
                             note="gene fixed effects; residue type as covariate; gene-cluster bootstrap; " + NOTE))
for d in ("D0", "D1"):
    rows.append(dict(block="G22", task="amend", signal="S4", signal_label=SIG["S4"], population="P-S/T (post-hoc)", design=d, index="V50",
                     tier_rule="study-cluster rule (884)", variation="post-hoc: serine and threonine records only; gene-cluster bootstrap", decision_bearing="no",
                     status="planned", note="restriction of the G1 population (P-all) to serine and threonine; " + NOTE))
n0 = len(m)
for i, r in enumerate(rows, 1):
    r["contrast_id"] = f"{r['block']}-{n0 + i:04d}"
with open(path, "a", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    for r in rows:
        w.writerow({c: r[c] for c in cols})
print(n0, "->", n0 + len(rows), sum(r["block"] == "G21" for r in rows), "G21 rows,", sum(r["block"] == "G22" for r in rows), "G22 rows")
