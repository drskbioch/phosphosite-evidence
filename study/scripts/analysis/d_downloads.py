"""Rev 6 downloads (Tasks D3, D4, E, H2). Files go to study/analysis/data_private/dl (git-ignored).

Records source URL, version stamp, retrieval date, size and SHA-256 for every file in
results/downloads_manifest.json. A failed or incomplete download is recorded with the exact error and
does not stop the others. Usage: python d_downloads.py [pharos goa paxdb atlas idmap]
"""
import datetime
import hashlib
import json
import sys

import pandas as pd
import requests

from common import PRIV, RES

DL = PRIV / "dl"
DL.mkdir(parents=True, exist_ok=True)
MAN = RES / "downloads_manifest.json"
man = json.load(open(MAN)) if MAN.exists() else {}
today = datetime.date.today().isoformat()


def save(key, rec):
    man[key] = rec
    json.dump(man, open(MAN, "w"), indent=1)
    print(key, json.dumps(rec)[:400], flush=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def fetch(url, dest, key, version=None, expect_md5=None):
    """Streaming download; size is compared with Content-Length."""
    rec = {"source_url": url, "retrieved": today, "version": version}
    try:
        with requests.get(url, stream=True, timeout=(30, 120)) as r:
            r.raise_for_status()
            rec["last_modified"] = r.headers.get("Last-Modified")
            want = int(r.headers.get("Content-Length", -1))
            tmp = dest.with_name(dest.name + ".part")
            with open(tmp, "wb") as f:
                for chunk in r.iter_content(1 << 20):
                    f.write(chunk)
        got = tmp.stat().st_size
        if want >= 0 and got != want:
            rec["error"] = f"incomplete: {got} of {want} bytes"
            save(key, rec)
            return rec
        tmp.replace(dest)
        rec.update(file=str(dest.relative_to(PRIV.parent.parent)), bytes=got, content_length=want, sha256=sha256(dest))
        if expect_md5:
            import hashlib as h
            m = h.md5(open(dest, "rb").read()).hexdigest()
            rec["md5"] = m
            rec["md5_expected"] = expect_md5
            if m != expect_md5:
                rec["error"] = f"md5 mismatch {m} != {expect_md5}"
    except Exception as e:  # noqa: BLE001
        rec["error"] = f"{type(e).__name__}: {e}"
    save(key, rec)
    return rec


def pharos():
    U = "https://pharos-api.ncats.io/graphql"

    def q(s):
        r = requests.post(U, json={"query": s}, timeout=180)
        r.raise_for_status()
        j = r.json()
        assert "errors" not in j, j["errors"]
        return j["data"]

    ver = q("{ dbVersion }")["dbVersion"]
    tot = q("{ targets { count } }")["targets"]["count"]
    rows, skip = [], 0
    while skip < tot:
        d = q("{ targets { targets(top:1000, skip:%d){ uniprot sym tdl fam } } }" % skip)["targets"]["targets"]
        if not d:
            break
        rows += d
        skip += len(d)
    df = pd.DataFrame(rows)
    rec = {"source_url": U, "version": ver, "retrieved": today, "count_reported": tot, "rows": len(df),
           "distinct_uniprot": int(df["uniprot"].nunique()), "tdl_counts": df["tdl"].value_counts().to_dict(),
           "query": "targets{targets(top,skip){uniprot sym tdl fam}}"}
    if len(df) != tot or df["uniprot"].nunique() != len(df.dropna(subset=["uniprot"]).drop_duplicates("uniprot")):
        rec["error"] = f"rows {len(df)} vs reported {tot}"
    df.to_parquet(PRIV / "pharos_tdl.parquet")
    save("pharos_tdl", rec)


def goa():
    base = "https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz"
    md5 = None
    try:
        txt = requests.get(base + ".md5", timeout=60).text.split()[0]
        md5 = txt if len(txt) == 32 and all(c in "0123456789abcdef" for c in txt) else None
    except Exception:  # noqa: BLE001
        pass
    fetch(base, DL / "goa_human.gaf.gz", "goa_human", version=None, expect_md5=md5)
    # version stamp from the file header (!gaf-version / !date-generated lines)
    import gzip
    hdr = []
    try:
        with gzip.open(DL / "goa_human.gaf.gz", "rt") as f:
            for line in f:
                if not line.startswith("!"):
                    break
                hdr.append(line.strip())
        man["goa_human"]["version"] = [h for h in hdr if "generated" in h.lower() or "version" in h.lower()][:4]
        save("goa_human", man["goa_human"])
    except Exception as e:  # noqa: BLE001
        man["goa_human"]["header_error"] = str(e)
        save("goa_human", man["goa_human"])


def paxdb():
    for t in (9606, 10090, 10116):
        url = f"https://pax-db.org/downloads/6.0/datasets/{t}/{t}-WHOLE_ORGANISM-integrated.txt"
        dest = DL / f"px_{t}.txt"
        rec = fetch(url, dest, f"paxdb_{t}", version="PaxDb 6.0 (v6.0 download path)")
        if "error" not in rec:
            hdr = [ln.strip() for ln in open(dest, encoding="utf-8", errors="replace") if ln.startswith("#")][:12]
            rec["header"] = hdr
            save(f"paxdb_{t}", rec)


def atlas():
    url = "https://peptideatlas.org/builds/human/phospho/202512/atlas_build_606.tsv.zip"
    fetch(url, DL / "atlas_build_606.tsv.zip", "peptideatlas_phospho_606", version="phospho build 606 (directory 202512)")


def idmap():
    try:
        rel = requests.get("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt", timeout=60).text.strip().splitlines()[:3]
    except Exception as e:  # noqa: BLE001
        rel = f"reldate.txt failed: {e}"
    for sp in ("HUMAN_9606", "MOUSE_10090", "RAT_10116"):
        url = f"https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/{sp}_idmapping.dat.gz"
        fetch(url, DL / f"{sp}_idmapping.dat.gz", f"uniprot_idmapping_{sp}", version=rel)


if __name__ == "__main__":
    todo = sys.argv[1:] or ["pharos", "goa", "paxdb", "atlas", "idmap"]
    for t in todo:
        try:
            globals()[t]()
        except Exception as e:  # noqa: BLE001
            save(t, {"error": f"{type(e).__name__}: {e}", "retrieved": today})
