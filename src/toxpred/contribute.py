"""Turn a file of molecules and measurements into a record others can pool.

Toxicity findings are the least sensitive data a discovery organisation holds.
A compound that shows toxicity is usually deprioritised or its program
terminated, so it is chemistry nobody will prosecute a claim on or advance, and
a finding that cost one company a program can stop three others repeating it.
That is the data a shared reference set most needs, and it is the data its owner
has the least reason to hold back. This writes the contribution in one of two
forms.

    fingerprint   the Morgan fingerprint and the endpoint values. Carries both
                  signals, structural neighbours and cross-family reach, and is
                  the better contribution once the chemistry is settled, which
                  for a toxicity finding it usually is.

    counts        four integers per compound, computed locally against the
                  public Reverse Screen index: protein families reached,
                  proteins reached, indexed ligands matched, and safety panel
                  proteins reached. Not a structural descriptor; no published
                  method recovers a structure from it. Carries the cross-family
                  signal only.

Nothing is uploaded. The output is a file you can inspect before you send it
anywhere.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

SCHEMA = 1
FORMATS = ("fingerprint", "counts")


def _read_input(path):
    """-> list of {smiles, endpoints{}}. Accepts CSV or TSV with a header.

    One column must be named smiles. Every other column is treated as an
    endpoint measurement, and blank cells are skipped rather than recorded as
    a value.
    """
    path = Path(path)
    text = path.read_text()
    sniff = "\t" if text.count("\t") > text.count(",") else ","
    rows = list(csv.DictReader(text.splitlines(), delimiter=sniff))
    if not rows:
        raise SystemExit("%s has no rows" % path)
    cols = {c.lower(): c for c in rows[0]}
    if "smiles" not in cols:
        raise SystemExit("%s needs a column named smiles; found: %s"
                         % (path, ", ".join(rows[0])))
    smi_col = cols["smiles"]
    out = []
    for r in rows:
        smi = (r.get(smi_col) or "").strip()
        if not smi:
            continue
        ep = {k: v.strip() for k, v in r.items()
              if k != smi_col and v is not None and v.strip() != ""}
        out.append(dict(smiles=smi, endpoints=ep))
    return out


def _counts(records, home, cutoff=0.5, family_map=None):
    """Screen locally and keep only the four integers."""
    from . import data as D
    from .reach import PharmCastMissing, fetch_model, read_family_map, reach

    files = D.fetch_reverse_screen(home)
    model = fetch_model(home)
    fam = read_family_map(family_map) if family_map else None
    rows = reach([r["smiles"] for r in records], files["index"], files["sites"],
                 model, family_map=fam, cutoff=cutoff)
    for rec, got in zip(records, rows):
        rec["counts"] = {
            "families_reached": got.get("families_reached"),
            "proteins_reached": got["targets_reached"],
            "matches": got["matches"],
        }
        rec.pop("smiles")
    return files["version"], model.name


def _fingerprints(records, radius=2, nbits=1024):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nbits)
    kept = []
    for rec in records:
        m = Chem.MolFromSmiles(rec["smiles"])
        if m is None:
            continue
        on = [int(b) for b in gen.GetFingerprint(m).GetOnBits()]
        kept.append(dict(bits=on, endpoints=rec["endpoints"]))
    return kept


def build(input_path, out_path, fmt="counts", home=None, cutoff=0.5,
          family_map=None, radius=2, nbits=1024):
    """Write a contribution file. Returns a short summary dict."""
    if fmt not in FORMATS:
        raise SystemExit("format must be one of: %s" % ", ".join(FORMATS))
    records = _read_input(input_path)
    meta = dict(schema=SCHEMA, format=fmt, n_compounds=len(records))

    if fmt == "counts":
        from .data import cache_dir
        version, model = _counts(records, cache_dir(home), cutoff, family_map)
        meta.update(index_version=version, fingerprint=model,
                    similarity_cutoff=cutoff,
                    note="four integers per compound. Not a structural "
                         "descriptor and no structure is included, so this "
                         "format is usable for chemistry still in play.")
        payload = records
    else:
        payload = _fingerprints(records, radius, nbits)
        meta.update(fingerprint="Morgan", radius=radius, bits=nbits,
                    n_compounds=len(payload),
                    note="carries both signals. Published work recovers a "
                         "fraction of structures from folded fingerprints, so "
                         "this format suits settled chemistry, which a toxicity "
                         "finding usually is.")

    endpoints = sorted({k for r in payload for k in r["endpoints"]})
    meta["endpoints"] = endpoints
    out = Path(out_path)
    out.write_text(json.dumps(dict(meta=meta, records=payload), indent=1))
    meta["written_to"] = str(out)
    return meta
