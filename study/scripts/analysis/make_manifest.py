"""Write study/analysis/grid_manifest.csv from section 3 of the analysis plan (one row per contrast)."""
import csv
from common import REV6

SIG = {"S1": "domain overlap", "S2": "conservation", "S3": "disorder", "S4": "kinase motif (sequence)", "S5": "kinase motif (lookup)"}
rows = []


def add(block, task, signal, pop, design, index="V50", tier="study-cluster rule (884)", extra="", bearing="yes", status="planned", note=""):
    rows.append(dict(block=block, task=task, signal=signal, signal_label=SIG[signal], population=pop, design=design,
                     index=index, tier_rule=tier, variation=extra, decision_bearing=bearing, status=status, note=note))


def bear(s):
    return "no" if s == "S5" else "yes"


for s in SIG:
    for d in ("D0", "D1", "D2"):
        add("G1", "B,C", s, "P-all" + (" (human)" if s == "S1" else ""), d, bearing=bear(s),
            note="database lookup: reported, not decision-bearing" if s == "S5" else "")
for s in ("S2", "S3", "S4", "S5"):
    for d in ("D0", "D1", "D2"):
        add("G2", "B,C", s, "P-human", d, bearing=bear(s))
for d in ("D0", "D1", "D2"):
    add("G3", "B,C", "S1", "P-annot", d, bearing="no", note="diagnostic only (D6): population defined from the outcome")
for s in ("S1", "S2", "S3", "S4"):
    for v in ("V0", "V20", "V100", "VT2025", "VT2026"):
        add("G4", "D", s, "P-all", "D3", index=v, extra="index variant", note="VT2026 is the reference-year sensitivity" if v == "VT2026" else "")
    for d in ("D0", "D1", "D5", "D4"):
        add("G5", "E", s, "P-abund (human)", d)
    for sp in ("mouse", "rat"):
        for d in ("D0", "D1", "D5", "D4"):
            add("G5", "E", s, f"P-abund ({sp})", d, status="conditional", note="run only if abundance data exist and strata can be filled; otherwise reported as not run")
    for r in ("S", "T", "Y"):
        for d in ("D0", "D1"):
            add("G6", "F", s, f"P-{r}", d)
    for k in [f"drop_cluster_{i}" for i in range(1, 11)] + ["drop_consortium_7_PMIDs"]:
        for d in ("D0", "D1"):
            add("G7", "F", s, "P-all", d, extra=k)
    for gd in (5, 10, 20):
        for ga in (2, 3):
            for d in ("D0", "D1"):
                add("G8", "F", s, "P-all", d, tier=f"GOLD depth>={gd} or DAS>={ga}; SILVER DAS>=1 or depth>=3", extra=f"gold_depth={gd};gold_das={ga}")
    for m in ("PSP LT-supported vs MS-only", "PSP evidence category high (n>=5) vs low (n=1)", "PeptideAtlas confidence categories", "UniProt MOD_RES experimental vs inferred ECO"):
        for d in ("D0", "D1", "D4"):
            cond = m.startswith(("PeptideAtlas", "UniProt"))
            add("G9", "H", s, "PSP∩site table (human)" if m.startswith("PSP") else "measure-specific", d, extra=m,
                status="conditional" if cond else "planned",
                note=("stop if not downloadable" if m.startswith("PeptideAtlas") else "optional") if cond else "mouse/rat only if intersection allows strata")
    add("G10", "I", s, "PSP human sites, proteins with both kinds", "within-protein fixed effects", index="none (within protein)", extra="LT-supported vs MS-only")
    add("G12", "B,C", s, "P-all", "D0", tier="depth-only (GOLD depth>=10, SILVER 3-9, BRONZE<=2)")
    add("G12", "B,C", s, "P-all", "D1", tier="depth-only (GOLD depth>=10, SILVER 3-9, BRONZE<=2)")
    for c in (500, 1000):
        for d in ("D0", "D1"):
            add("G13", "B,C", s, "P-all", d, tier=f"study-cluster rule ({c})", extra=f"cluster_threshold={c}",
                note="cluster merges and tiers identical to 884 in b_cluster_rule_info.json")
for s in ("S1", "S2", "S3", "S4", "S5"):
    for d in ("R5_CEM_primary_index", "R5_CEM_secondary_index", "R5_propensity"):
        add("G14", "B", s, "as in Rev 5", d, index="V0 (original)", tier="set-A tier (S1); set-B tier (others)", bearing="no", extra="protein-cluster CI")
for d in ("D0", "D1", "D2"):
    add("G15", "B", "S1", "human", d, tier="set-A tier", bearing="no", note="labelled comparison row for D1")
G11 = [("DeepPhos", "S/T", 2019), ("DeepPhos", "S/T", 2018), ("MusiteDeep", "S/T", 2020), ("MusiteDeep", "S/T", 2017), ("PhosF3C", "Y", 2024)]
for tool, res, yr in G11:
    rows.append(dict(block="G11", task="J", signal="AUROC gap", signal_label="top-three minus bottom-three decile AUROC gap", population=f"{tool} {res} sites, positives first reported after {yr}",
                     design="equivalence test (margin 0.05), protein-cluster bootstrap", index="decile of study effort", tier_rule="n/a", variation=f"restriction_year={yr}",
                     decision_bearing="yes", status="planned", note="restriction years from earlier_results/tool_versions.json; PhosF3C tyrosine only; no pooled ranking"))
    rows.append(dict(block="G11", task="J", signal="AUROC gap", signal_label="top-three minus bottom-three decile AUROC gap", population=f"{tool} {res} sites, unrestricted",
                     design="equivalence test (margin 0.05), protein-cluster bootstrap", index="decile of study effort", tier_rule="n/a", variation="unrestricted",
                     decision_bearing="yes", status="planned", note=""))
seen = set(); out = []
for r in rows:
    key = tuple(r.values())
    if key in seen:
        continue
    seen.add(key); out.append(r)
for i, r in enumerate(out, 1):
    r["contrast_id"] = f"{r['block']}-{i:04d}"
cols = ["contrast_id"] + [c for c in out[0] if c != "contrast_id"]
with open(REV6 / "grid_manifest.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(out)
from collections import Counter
print(len(out), dict(Counter(r["block"] for r in out)))
