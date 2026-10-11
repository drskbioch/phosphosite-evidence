"""Write the title of the work into README.md, CITATION.cff and .zenodo.json from the single file TITLE.txt.

Usage: python scripts/set_title.py            rewrite the three files
       python scripts/set_title.py --check    change nothing; exit status 1 if a file does not carry the title of TITLE.txt
README.md: every text between <!--title--> and <!--/title--> is replaced. CITATION.cff: the line `title:`. .zenodo.json: the key `title`, which is the title
with the prefix below. Only the standard library is used.
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ZENODO_PREFIX = "Analysis code and results for: "


def new_texts(title):
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    readme = re.sub(r"(<!--title-->).*?(<!--/title-->)", lambda m: m.group(1) + title + m.group(2), readme, flags=re.S)
    cff = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    quoted = title.replace(chr(92), chr(92) * 2).replace(chr(34), chr(92) + chr(34))
    cff = re.sub(r"(?m)^title:.*$", lambda m: "title: " + chr(34) + quoted + chr(34), cff, count=1)
    zen = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))
    zen["title"] = ZENODO_PREFIX + title
    return {"README.md": readme, "CITATION.cff": cff, ".zenodo.json": json.dumps(zen, indent=2, ensure_ascii=False) + chr(10)}


def main():
    title = (ROOT / "TITLE.txt").read_text(encoding="utf-8").strip()
    if not title or chr(10) in title:
        sys.exit("TITLE.txt must hold one non-empty line")
    check = "--check" in sys.argv
    changed = []
    for name, text in new_texts(title).items():
        p = ROOT / name
        if p.read_text(encoding="utf-8") != text:
            changed.append(name)
            if not check:
                p.write_text(text, encoding="utf-8", newline=chr(10))
    print(("would change: " if check else "changed: ") + (", ".join(changed) if changed else "nothing (all three files carry the title)"))
    sys.exit(1 if changed and check else 0)


if __name__ == "__main__":
    main()
