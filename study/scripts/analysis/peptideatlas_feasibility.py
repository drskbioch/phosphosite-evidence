"""Feasibility check, PeptideAtlas build 606: can site-level phosphosite evidence be reconstructed from the peptide tables?
Descriptive only; no evidence category is derived. Reads data_private/dl/atlas_build_606.tsv.zip and results/peptideatlas_606_inspection.json.
Output: results/peptideatlas_606_feasibility.json
"""
import collections
import json
import re
import zipfile

from common import PRIV, RES

NULL = chr(92) + "N"
Z = zipfile.ZipFile(PRIV / "dl" / "atlas_build_606.tsv.zip")
insp = json.load(open(RES / "peptideatlas_606_inspection.json"))

# 1. columns whose names suggest a localisation score or site-level confidence, over all tables
pat = re.compile(r"local|ptm|prophet|ascore|flr|site|conf|categ|prob|score", re.I)
hits = {t.replace(".tsv", ""): [c for c in cols if pat.search(c)] for t, cols in insp["tables"].items()}
hits = {t: c for t, c in hits.items() if c}

# 2. modified_peptide_instance: modification tokens, residues carrying [Phospho], placement ambiguity, null shares
R = re.compile(r"([A-Z])((?:\[[^\]]+\])*)")
res_count = collections.Counter()
place = collections.defaultdict(set)
n = n_ph = 0
nulls = collections.Counter()
cols_check = ["n_observations", "n_adjusted_observations", "best_probability", "best_adjusted_probability", "sample_ids", "search_batch_ids"]
with Z.open("modified_peptide_instance.tsv") as f:
    hdr = next(f).decode().rstrip("\n").split("\t")
    ix = {h: i for i, h in enumerate(hdr)}
    for line in f:
        p = line.decode("utf8", "replace").rstrip("\n").split("\t")
        n += 1
        for c in cols_check:
            if p[ix[c]] == NULL:
                nulls[c] += 1
        s = p[ix["modified_peptide_sequence"]]
        if "[Phospho]" not in s:
            continue
        n_ph += 1
        s2 = re.sub(r"^\[[^\]]+\]-", "", s)
        pos, strip = [], []
        for i, (r, m) in enumerate(R.findall(s2)):
            strip.append(r)
            if "[Phospho]" in m:
                res_count[r] += 1
                pos.append(i)
        place[(p[ix["peptide_instance_id"]], "".join(strip), len(pos))].add(tuple(pos))
amb = sum(1 for v in place.values() if len(v) > 1)
tot_tok = sum(res_count.values())

# 3. peptide_instance: protein-mapping multiplicity
nmap = collections.Counter()
with Z.open("peptide_instance.tsv") as f:
    hdr = next(f).decode().rstrip("\n").split("\t")
    ix = {h: i for i, h in enumerate(hdr)}
    for line in f:
        p = line.decode("utf8", "replace").rstrip("\n").split("\t")
        v = p[ix["n_protein_mappings"]]
        nmap["null" if v == NULL else ("1" if v == "1" else "2 or more")] += 1

out = {
    "columns_matching_localisation_or_site_confidence_pattern": hits,
    "modified_peptide_instances": n,
    "with_phospho_token": n_ph,
    "null_counts_in_modified_peptide_instance": dict(nulls),
    "phospho_tokens_by_residue": dict(res_count.most_common()),
    "phospho_tokens_total": tot_tok,
    "share_on_residues_other_than_STY": sum(v for k, v in res_count.items() if k not in "STY") / tot_tok,
    "peptide_by_phospho_count_groups": len(place),
    "groups_with_more_than_one_placement": amb,
    "share_groups_with_more_than_one_placement": amb / len(place),
    "peptide_instance_n_protein_mappings": dict(nmap),
    "build": {"is_ptm_build": "Y", "description": "Human Phosphoproteome 2025-12; phosphopeptide enrichment samples (atlas_build.tsv)"},
}
json.dump(out, open(RES / "peptideatlas_606_feasibility.json", "w"), indent=1)
print(json.dumps(out, indent=1))
