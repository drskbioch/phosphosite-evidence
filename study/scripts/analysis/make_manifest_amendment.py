"""Append the blocks added by the 2026-10-01 plan amendment to grid_manifest.csv (existing rows are not changed).
Run once, before any contrast of these blocks is run. Idempotent: does nothing if G16 rows already exist."""
import pandas as pd
import csv
from common import REV6

path = REV6 / "grid_manifest.csv"
m = pd.read_csv(path).fillna("")
if (m["block"] == "G16").any():
    raise SystemExit("amendment rows already present")
cols = list(m.columns)
SIG = {"S1": "domain overlap", "S2": "conservation", "S3": "disorder", "S4": "kinase motif (sequence)", "S5": "kinase motif (lookup)"}
AMEND = "added by the plan amendment of 2026-10-01 (after the interim report)"
rows = []


def add(block, task, s, pop, design, tier="study-cluster rule (884)", index="V50", var="", bearing="yes", note=""):
    rows.append(dict(block=block, task=task, signal=s, signal_label=SIG[s], population=pop, design=design, index=index, tier_rule=tier,
                     variation=var, decision_bearing=bearing, status="planned", note=(note + "; " if note else "") + AMEND))


def bear(s):
    return "no" if s == "S5" else "yes"


# G16 canonical accessions (G1 and G2 populations)
for s in SIG:
    for d in ("D0", "D1"):
        add("G16", "amend", s, "P-all (human)" if s == "S1" else "P-all, Swiss-Prot reviewed canonical", d, var="canonical accessions only", bearing=bear(s))
for s in ("S2", "S3", "S4", "S5"):
    for d in ("D0", "D1"):
        add("G16", "amend", s, "P-human, Swiss-Prot reviewed canonical", d, var="canonical accessions only", bearing=bear(s))
# fix the S1 population label
for r in rows:
    if r["signal"] == "S1":
        r["population"] = "P-all (human), Swiss-Prot reviewed canonical"
# G17 gene-cluster bootstrap of G1, G2, G3, G5 (human) D0 and D1
src = m[m["block"].isin(["G1", "G2", "G3"]) & m["design"].isin(["D0", "D1"])]
for _, r in src.iterrows():
    add("G17", "amend", r["signal"], r["population"], r["design"], bearing=r["decision_bearing"], var=f"gene-cluster bootstrap; source block {r['block']}",
        note=f"same records as {r['contrast_id']}")
src5 = m[(m["block"] == "G5") & m["design"].isin(["D0", "D1"]) & (m["population"] == "P-abund (human)")]
for _, r in src5.iterrows():
    add("G17", "amend", r["signal"], r["population"], r["design"], bearing=r["decision_bearing"], var="gene-cluster bootstrap; source block G5",
        note=f"same records as {r['contrast_id']}")
# G18 post-hoc matcher D2b: same 11 contrasts as the D2 rows of G1, G2, G3, G15
src = m[m["block"].isin(["G1", "G2", "G3", "G15"]) & (m["design"] == "D2")]
for _, r in src.iterrows():
    add("G18", "amend", r["signal"], r["population"], "D2b", tier=r["tier_rule"], bearing=r["decision_bearing"],
        var=f"post-hoc matcher: 1:1 without replacement, caliper 0.2 SD, S/T/Y only, random order, >= 20 order replicates; source block {r['block']}",
        note=f"counterpart of {r['contrast_id']}")
# G19 tier-definition sensitivity
TIERS = ["set-B tiers (primary)", "set-A tiers", "depth-only tiers", "agreement-only (agreement >= 2 vs 0)"]
for s in ("S1", "S2", "S3", "S4"):
    for t in TIERS:
        for d in ("D0", "D1"):
            add("G19", "amend", s, "P-human", d, tier=t, var="tier-definition sensitivity")
for s in ("S2", "S3", "S4"):
    for t in (TIERS[0], TIERS[2], TIERS[3]):
        for d in ("D0", "D1"):
            add("G19", "amend", s, "P-all", d, tier=t, var="tier-definition sensitivity", note="set-A tiers exist for human sites only")
# G20 diagnostics
add("G20", "amend", "S1", "all D1 and D3 rows run", "pre-matching balance", var="SMD of the primary index, GOLD vs BRONZE, before matching", bearing="no",
    note="one row; the value is attached to every D1 and D3 contrast of the results tables")
for s in ("S1", "S2", "S3", "S4", "S5"):
    for pop in (["P-all (human)"] if s == "S1" else ["P-all", "P-human"]):
        for tier in ("GOLD", "BRONZE"):
            add("G20", "amend", s, pop, "index-outcome association", var=f"within {tier} records; gene-cluster interval", bearing="no")
add("G20", "amend", "S1", "protein level", "index vs sites per protein", var="Spearman, gene-cluster interval", bearing="no")
add("G20", "amend", "S1", "protein level", "index vs LT_LIT-supported sites per gene", var="Spearman, gene-cluster interval", bearing="no", note="requires the PhosphoSitePlus file")
for t in ("annotated", "not annotated"):
    add("G20", "amend", "S1", f"P-all (human), proteins with {'at least one' if t == 'annotated' else 'no'} domain-annotated site", "domain-annotation coverage strata",
        var="crude odds ratio and tier shares per stratum", bearing="no", note="in the stratum without an annotated site the outcome is constant, so the odds ratio is undefined")
n0 = len(m)
out = pd.DataFrame(rows)
for b, g in out.groupby("block", sort=False):
    pass
ids = []
cnt = {}
for r in rows:
    cnt[r["block"]] = cnt.get(r["block"], 0) + 1
    ids.append(f"{r['block']}-{n0 + len(ids) + 1:04d}")
for r, i in zip(rows, ids):
    r["contrast_id"] = i
add_df = pd.DataFrame(rows)[cols]
with open(path, "a", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    for _, r in add_df.iterrows():
        w.writerow(r.to_dict())
print(len(m), "->", len(m) + len(add_df), cnt)
