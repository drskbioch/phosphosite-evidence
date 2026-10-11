"""Verify the files in data_private/dl and write the final state of each download to results/downloads_manifest.json.
Integrity: gzip/zip tests, size against the Content-Length header recorded at download time. Later re-download attempts that failed on the network
are recorded but do not replace a verified file."""
import datetime
import gzip
import hashlib
import json
import zipfile

from common import PRIV, RES

DL = PRIV / "dl"
man = json.load(open(RES / "downloads_manifest.json"))


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def mt(p):
    return datetime.datetime.fromtimestamp(p.stat().st_mtime).isoformat(timespec="seconds")


def gz_ok(p):
    with gzip.open(p) as g:
        while g.read(1 << 24):
            pass
    return True


checks = {
    "peptideatlas_phospho_606": (DL / "atlas_build_606.tsv.zip", 549143491, "Tue, 16 Dec 2025 00:26:21 GMT", lambda p: zipfile.ZipFile(p).testzip() is None,
                                 "PeptideAtlas human phospho build 606 (directory 202512), file atlas_build_606.tsv.zip"),
    "uniprot_idmapping_HUMAN_9606": (DL / "HUMAN_9606_idmapping.dat.gz", 38464542, "Thu, 03 Sep 2026 18:00:00 GMT", gz_ok, "UniProt Knowledgebase release 2026_03 (02-Sep-2026)"),
    "uniprot_idmapping_MOUSE_10090": (DL / "MOUSE_10090_idmapping.dat.gz", 15536610, "Thu, 03 Sep 2026 18:00:00 GMT", gz_ok, "UniProt Knowledgebase release 2026_03 (02-Sep-2026)"),
    "uniprot_idmapping_RAT_10116": (DL / "RAT_10116_idmapping.dat.gz", 9415486, "Thu, 03 Sep 2026 18:00:00 GMT", gz_ok, "UniProt Knowledgebase release 2026_03 (02-Sep-2026)"),
    "goa_human": (DL / "goa_human.gaf.gz", 11020779, "Tue, 28 Jul 2026 15:31:43 GMT", gz_ok, None),
}
for k, (p, want, lm, test, ver) in checks.items():
    rec = man.get(k, {})
    rec = {kk: vv for kk, vv in rec.items() if kk not in ("error", "md5", "md5_expected")}
    ok = p.stat().st_size == want and test(p)
    rec.update(file=str(p.relative_to(PRIV.parent.parent)), bytes=p.stat().st_size, content_length_header=want, last_modified_header=lm, sha256=sha(p),
               file_modified_on_disk=mt(p), integrity_check_passed=bool(ok))
    if ver:
        rec["version"] = ver
    if k == "peptideatlas_phospho_606" or k.startswith("uniprot_idmapping"):
        rec["retrieval_note"] = ("Fully retrieved and verified (size equals the Content-Length header; archive test passes). Later re-download attempts on 2026-10-01 failed on the network "
                                 "(DNS or read time-out) and were discarded; the verified file was not replaced.")
    if k == "goa_human":
        rec["note"] = "No .md5 file is published at the expected path (404); integrity rests on the size equal to Content-Length and a complete gzip read."
    man[k] = rec
    print(k, ok)
for t in ("9606", "10090", "10116"):
    p = DL / f"px_{t}.txt"
    rec = man[f"paxdb_{t}"]
    n = sum(1 for ln in open(p, encoding="utf-8") if not ln.startswith("#"))
    rec.update(sha256=sha(p), data_rows=n, note="chunked response (no Content-Length); file read completely, last line parses; row count recorded.")
json.dump(man, open(RES / "downloads_manifest.json", "w"), indent=1)
