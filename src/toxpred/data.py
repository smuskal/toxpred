"""Fetch the public reference data. Nothing here is proprietary.

Three sources, all freely downloadable:

    CATMoS   rat oral acute toxicity, curated by NICEATM and US EPA and
             distributed with the OPERA toolset. Experimental values only are
             kept; the same file also carries consensus model predictions and
             those are never used as reference values.
    TOXRIC   thirty further standardized toxicity endpoints.
    index    the Reverse Screen index of ligands solved in protein structures,
             used only for the optional cross-family reach signal.

Everything lands under a cache directory you control, by default
~/.toxpred, so nothing is written into your environment or your project.
"""
from __future__ import annotations

import csv
import os
import re
import zipfile
from pathlib import Path

import requests

OPERA_DATA = ("https://raw.githubusercontent.com/NIEHS/OPERA/master/"
              "OPERA_Data.zip")
CATMOS_MEMBER = "CATMoS_QR50k.sdf"
TOXRIC_30 = "https://ndownloader.figshare.com/files/49694949"
REVERSE_SCREEN_API = "https://reversescreen.ai/api/download"

CATMOS_FIELDS = ["DTXSID", "casrn", "Name", "Canonical_QSARr",
                 "InChI Key_QSARr", "CATMoS_logLD50_data", "CATMoS_LD50_mgkg",
                 "CATMoS_EPA_data", "CATMoS_GHS_data", "CATMoS_VT_data",
                 "CATMoS_NT_data", "CATMoS_Tr_Tst"]
CATMOS_COLS = ["dtxsid", "casrn", "name", "smiles", "inchikey", "log_ld50",
               "ld50_mgkg", "epa_cat", "ghs_cat", "very_toxic", "non_toxic",
               "split"]


def cache_dir(path: str | os.PathLike | None = None) -> Path:
    d = Path(path or os.environ.get("TOXPRED_HOME", Path.home() / ".toxpred"))
    d.mkdir(parents=True, exist_ok=True)
    return d


def _download(url: str, dest: Path, desc: str = "") -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    import sys as _sys
    show = _sys.stdout.isatty()
    with requests.get(url, stream=True, timeout=600) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        done = 0
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
                done += len(chunk)
                if total and show:
                    print("\r  %s %5.1f%%" % (desc or dest.name,
                                              100.0 * done / total),
                          end="", flush=True)
    print(("\r" if show else "  ") + "%s done" % (desc or dest.name))
    tmp.rename(dest)
    return dest


def _sdf_records(path: Path):
    cur, field = {}, None
    with open(path, errors="ignore") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("$$$$"):
                yield cur
                cur, field = {}, None
            elif line.startswith("> <") or line.startswith(">  <"):
                m = re.match(r"> *<([^>]*)>", line)
                field = m.group(1) if m else None
            elif field is not None:
                if line.strip():
                    cur[field] = line.strip()
                field = None


def fetch_catmos(home: Path) -> Path:
    """-> a csv of the EXPERIMENTAL rat oral acute toxicity records."""
    out = home / "catmos_experimental.csv"
    if out.exists():
        return out
    zpath = _download(OPERA_DATA, home / "OPERA_Data.zip", "OPERA data")
    sdf = home / CATMOS_MEMBER
    if not sdf.exists():
        with zipfile.ZipFile(zpath) as z:
            z.extract(CATMOS_MEMBER, home)
    n = kept = 0
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(CATMOS_COLS)
        for rec in _sdf_records(sdf):
            n += 1
            # a record with no train/test flag was never measured; it is in the
            # file to carry a model prediction, and is not a reference value
            if not rec.get("CATMoS_Tr_Tst") or not rec.get("Canonical_QSARr"):
                continue
            w.writerow([rec.get(k, "") for k in CATMOS_FIELDS])
            kept += 1
    print("  CATMoS: %d records read, %d experimental rows kept" % (n, kept))
    return out


def fetch_toxric(home: Path) -> Path:
    """-> the directory holding the thirty endpoint tables."""
    out = home / "toxric"
    if out.exists() and any(out.glob("**/*.csv")):
        return out
    z = _download(TOXRIC_30, home / "toxric_30_datasets.zip", "TOXRIC")
    with zipfile.ZipFile(z) as zf:
        zf.extractall(out)
    return out


def reverse_screen_manifest() -> dict:
    r = requests.get(REVERSE_SCREEN_API, timeout=120)
    r.raise_for_status()
    return r.json()


def fetch_reverse_screen(home: Path) -> dict:
    """The published index, its ligand and site tables, checksum verified."""
    import hashlib
    man = reverse_screen_manifest()
    want = {}
    for f in man["files"]:
        if any(k in f["name"] for k in ("index_", "sites_", "targets_")):
            want[f["name"]] = f["sha256"]
    got = {}
    for name, sha in want.items():
        p = _download("%s/%s" % (REVERSE_SCREEN_API, name), home / name, name)
        h = hashlib.sha256()
        with open(p, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        if h.hexdigest() != sha:
            raise RuntimeError("checksum mismatch for %s" % name)
        kind = ("index" if "index_" in name
                else "sites" if "sites_" in name else "targets")
        got[kind] = p
    got["version"] = man["version"]
    print("  reverse screen index version %s, checksums verified"
          % man["version"])
    return got
