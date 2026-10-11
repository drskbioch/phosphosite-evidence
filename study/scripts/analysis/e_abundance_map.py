"""Task E1: PaxDb integrated whole-organism abundance (ppm) mapped to UniProt accessions.

PaxDb rows: gene name, STRING external id (<taxid>.ENSP...), abundance. Mapping: UniProt idmapping by_organism file, id type STRING.
If several PaxDb entries map to one accession, the maximum abundance is kept (rule fixed before any outcome was analysed; the number of
accessions with more than one entry is reported). Unmapped PaxDb entries and coverage of proteins and records are counted.
Output: data_private/abundance.parquet (uniprot_id, species, abundance_ppm) and results/e_abundance_coverage.json.
"""
import gzip
import json

import numpy as np
import pandas as pd

import base as B
from common import PRIV, RES

DL = PRIV / "dl"
b = B.load(["uniprot_id", "species", "SEI_V50", "tier", "mapped_to_gene"])
SP = {9606: ("human", "HUMAN_9606"), 10090: ("mouse", "MOUSE_10090"), 10116: ("rat", "RAT_10116")}
frames, cov = [], {}
for tax, (sp, idm) in SP.items():
    px = pd.read_csv(DL / f"px_{tax}.txt", sep="\t", comment="#", header=None, names=["gene_name", "string_id", "abundance_ppm"])
    px["abundance_ppm"] = pd.to_numeric(px["abundance_ppm"], errors="raise")
    assert (px["abundance_ppm"] > 0).all()
    mp = []
    with gzip.open(DL / f"{idm}_idmapping.dat.gz", "rt") as f:
        for line in f:
            a, t, v = line.rstrip("\n").split("\t")
            if t == "STRING":
                mp.append((a, v))
    mp = pd.DataFrame(mp, columns=["uniprot_id", "string_id"])
    m = px.merge(mp, on="string_id", how="left")
    n_unmapped = int(m.loc[~m["string_id"].isin(mp["string_id"]), "string_id"].nunique())
    m = m.dropna(subset=["uniprot_id"])
    dup = int((m.groupby("uniprot_id").size() > 1).sum())
    g = m.groupby("uniprot_id", as_index=False)["abundance_ppm"].max()
    g["species"] = sp
    frames.append(g)
    sub = b[b["species"] == sp]
    prot = sub.drop_duplicates("uniprot_id")
    has = prot["uniprot_id"].isin(g["uniprot_id"])
    rec_has = sub["uniprot_id"].isin(g["uniprot_id"])
    idx_ok = sub["SEI_V50"].notna()
    cov[sp] = {"paxdb_entries": int(len(px)), "paxdb_entries_unmapped_to_uniprot": n_unmapped, "uniprot_accessions_mapped": int(g["uniprot_id"].nunique()),
               "accessions_with_more_than_one_paxdb_entry": dup, "proteins_in_site_table": int(len(prot)), "proteins_with_abundance": int(has.sum()),
               "share_of_proteins": float(has.mean()), "records_in_site_table": int(len(sub)), "records_with_abundance": int(rec_has.sum()),
               "share_of_records": float(rec_has.mean()), "records_with_abundance_and_V50_index": int((rec_has & idx_ok).sum()),
               "proteins_with_abundance_and_V50_index": int(prot.loc[has & prot['uniprot_id'].isin(sub.loc[idx_ok, 'uniprot_id']), 'uniprot_id'].nunique()),
               "gold_bronze_records_with_abundance_and_index": int((rec_has & idx_ok & sub["tier"].isin(["GOLD", "BRONZE"])).sum())}
    print(sp, cov[sp], flush=True)
pd.concat(frames).to_parquet(PRIV / "abundance.parquet")
json.dump(cov, open(RES / "e_abundance_coverage.json", "w"), indent=1)
