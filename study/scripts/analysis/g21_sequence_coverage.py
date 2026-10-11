"""G21: share of records with a UniProt reviewed protein sequence (and matching residue), and the detectability indicator, by tier; per signal-by-population group.
Descriptive; output results/g21_sequence_coverage.json"""
import json

import g21_data as G
from common import RES

out = {}
for popl in ("P-all (S/T/Y)", "P-human (S/T/Y)"):
    for s in ("S2", "S3", "S4", "S5"):
        d, info = G.prepare(s, popl)
        r = {"records": int(len(d))}
        for tier, nm in ((1, "GOLD"), (0, "BRONZE")):
            x = d[d["gold"] == tier]
            h = x[x["has_seq"]]
            r[nm] = {"records": int(len(x)), "with_sequence_and_match": int(len(h)), "share_with_sequence": float(len(h) / len(x)), "share_detectable_among_with_sequence": float(h["det"].mean())}
        out[f"{s}|{popl}"] = r
        print(s, popl, r["GOLD"]["share_with_sequence"], r["BRONZE"]["share_with_sequence"], flush=True)
json.dump(out, open(RES / "g21_sequence_coverage.json", "w"), indent=1)
