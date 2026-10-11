"""Provenance of the numbers a report builder prints.

`jload(path)` reads a JSON result file and returns the same structure with every number replaced by a float or int subclass that remembers the file
and the key path it came from. A tracked number that is formatted (f-string, format(), or `fmt_int`) writes (text, file, key path, scale) to `LOG`.
Multiplying by a plain constant (100 * share) keeps the provenance with the scale; any other arithmetic returns a plain number, which is not tracked
(the checker then reports such a number as derived or as not bound). `RowList` attaches the log entries made while a table row was built to the row's
ID, and `dump_bindings` writes them to a sidecar file that `check_bindings.py` re-verifies against the result files.
"""
import json
from pathlib import Path

LOG = []
ROOT = Path(__file__).resolve().parents[2]  # study/


class T(float):
    """Tracked float."""

    def __new__(cls, v, src, scale=1.0):
        o = super().__new__(cls, v)
        o.src, o.scale = src, scale
        return o

    def __format__(self, spec):
        s = float.__format__(self, spec)
        LOG.append((s, self.src[0], self.src[1], self.scale))
        return s

    def __mul__(self, k):
        return T(float(self) * k, self.src, self.scale * k) if isinstance(k, (int, float)) and not isinstance(k, (T, TI, bool)) else float(self) * float(k)

    __rmul__ = __mul__

    def __neg__(self):
        return T(-float(self), self.src, -self.scale)

    def __abs__(self):
        return T(abs(float(self)), self.src, abs(self.scale))


class TI(int):
    """Tracked int."""

    def __new__(cls, v, src, scale=1):
        o = super().__new__(cls, v)
        o.src, o.scale = src, scale
        return o

    def __format__(self, spec):
        s = int.__format__(self, spec)
        LOG.append((s, self.src[0], self.src[1], self.scale))
        return s

    def __mul__(self, k):
        return T(int(self) * k, self.src, self.scale * k) if isinstance(k, (int, float)) and not isinstance(k, (T, TI, bool)) else int(self) * (float(k) if isinstance(k, float) else int(k))

    __rmul__ = __mul__

    def __neg__(self):
        return TI(-int(self), self.src, -self.scale)

    def __abs__(self):
        return TI(abs(int(self)), self.src, abs(self.scale))


def fmt_int(v):
    """f"{int(round(v)):,}" that keeps the provenance of a tracked number."""
    s = f"{int(round(v)):,}"
    if isinstance(v, (T, TI)):
        LOG.append((s, v.src[0], v.src[1], v.scale))
    return s


def track(o, rel, path=()):
    if isinstance(o, bool) or o is None:
        return o
    if isinstance(o, float):
        return T(o, (rel, ".".join(path)))
    if isinstance(o, int):
        return TI(o, (rel, ".".join(path)))
    if isinstance(o, dict):
        return {k: track(v, rel, path + (str(k),)) for k, v in o.items()}
    if isinstance(o, list):
        return [track(v, rel, path + (str(i),)) for i, v in enumerate(o)]
    return o


def jload(p):
    p = Path(p)
    try:
        rel = str(p.resolve().relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        rel = str(p)
    return track(json.load(open(p, encoding="utf-8")), rel)


class RowRec(list):
    def __init__(self, r, owner):
        super().__init__(r)
        self.owner = owner

    def __setitem__(self, i, v):
        super().__setitem__(i, v)
        self.owner.attach(self[0])


class RowList(list):
    """Rows of a table; the log entries made since the previous row (or the previous item assignment) belong to the row appended."""
    BIND = {}
    MARK = [0]

    def attach(self, rid):
        RowList.BIND.setdefault(str(rid), []).extend(LOG[RowList.MARK[0]:])
        RowList.MARK[0] = len(LOG)

    def append(self, r):
        self.attach(r[0])
        super().append(RowRec(r, self))


def dump_bindings(path):
    out = {rid: [list(e) for e in es] for rid, es in RowList.BIND.items()}
    Path(path).write_text(json.dumps(out, indent=0), encoding="utf-8")
