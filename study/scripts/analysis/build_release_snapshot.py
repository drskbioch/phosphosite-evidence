"""Build the archive snapshot: a zip of exactly the files listed in section 6 of RELEASE_SCOPE.md, and MANIFEST.sha256 (SHA-256 of every file).

Writes to study/analysis/data_private/release_snapshot/ (git-ignored, not part of any commit) unless --out DIR is given:
  release_snapshot.zip       the listed files, at their repository-relative paths, in sorted order, with fixed time stamps (the same files give the same zip)
  MANIFEST.sha256            one line per file: SHA-256, two spaces, path (readable by `sha256sum -c`)
  release_snapshot.zip.sha256  SHA-256 of the zip
  SNAPSHOT_COMMIT.txt        the commit that was checked out when the zip was built
The builder stops if a listed file is missing, has an excluded prefix (see check_release_scope.py), or is modified or untracked in git (the zip then equals the committed files).
Nothing is uploaded or deposited.
Usage: python build_release_snapshot.py [--out DIR]
"""
import hashlib
import re
import subprocess
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "study" / "scripts" / "release"))
from check_release_scope import EXCLUDED, listed_paths, nested_repo  # noqa: E402

OUT = Path(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv else ROOT / "study" / "analysis" / "data_private" / "release_snapshot"
FIXED_TIME = (2026, 1, 1, 0, 0, 0)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def git(*args):
    return subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True, encoding="utf-8").stdout


def main():
    rels = sorted(listed_paths())
    assert len(rels) == len(set(rels)), "duplicate path in the list"
    problems = []
    for r in rels:
        if not (ROOT / r).is_file():
            problems.append(f"missing: {r}")
        if r.startswith(tuple(p for p, _ in EXCLUDED)) or nested_repo(r):
            problems.append(f"excluded or in a nested repository: {r}")
    dirty = set()
    for line in git("status", "--porcelain", "--untracked-files=all", "--ignored=no").splitlines():
        dirty.add(line[3:].strip().strip('"'))
    tracked = set(git("ls-files").splitlines())
    for r in rels:
        if r in dirty:
            problems.append(f"modified or untracked in git: {r}")
        elif r not in tracked:
            problems.append(f"not tracked in git: {r}")
    if problems:
        print("\n".join(problems[:40]))
        print(f"{len(problems)} problem(s); no snapshot written")
        sys.exit(1)
    OUT.mkdir(parents=True, exist_ok=True)
    zpath = OUT / "release_snapshot.zip"
    lines = []
    with zipfile.ZipFile(zpath, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for r in rels:
            p = ROOT / r
            info = zipfile.ZipInfo(r, date_time=FIXED_TIME)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            z.writestr(info, p.read_bytes(), compresslevel=9)
            lines.append(f"{sha256(p)}  {r}")
    (OUT / "MANIFEST.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    zh = sha256(zpath)
    (OUT / "release_snapshot.zip.sha256").write_text(f"{zh}  release_snapshot.zip\n", encoding="utf-8", newline="\n")
    (OUT / "SNAPSHOT_COMMIT.txt").write_text(git("rev-parse", "HEAD"), encoding="utf-8", newline="\n")
    print(f"{len(rels)} files; zip {zpath} {zpath.stat().st_size:,} bytes; SHA-256 {zh}")
    print(f"commit {git('rev-parse', 'HEAD').strip()}")


if __name__ == "__main__":
    main()
