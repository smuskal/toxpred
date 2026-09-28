"""Cross-family reach: how many protein families a compound's chemistry reaches.

Every compound is compared to the ligands solved in protein structures, and the
proteins those matches were solved against are collected. Reach is the number of
distinct protein families among them. Compounds that reach widely are more often
toxic, and that is measured with no access to how many assay panels a compound
has been through, which is the confound that makes measured promiscuity hard to
read.

THIS FEATURE NEEDS PharmCast. The published index stores a pharmacophore
fingerprint per ligand, and a query has to be fingerprinted the same way to be
comparable. PharmCast is a separate open package:

    https://pharmcast.ai
    https://github.com/smuskal/pharmcast

Install it with `pip install git+https://github.com/smuskal/pharmcast.git`, or
skip reach entirely. Everything else in toxpred works without it.

Every fingerprint operation here is PharmCast's own: its reader for the index,
its batch predictor for the queries, and its similarity routine for the
comparison. None of the packing or bit ordering is reimplemented, so the two
sides cannot drift apart.
"""
from __future__ import annotations

import collections
from pathlib import Path

# The published PharmCast checkpoint. This is what the package downloads and
# what every number in this repository was produced with. Change it here only.
MODEL = "pharmcast_scp_v10.pt"
MODELS_BASE = "https://pharmcast.ai/models"
MODEL_URL = "%s/%s" % (MODELS_BASE, MODEL)
SUMS_URL = "%s/SHA256SUMS" % MODELS_BASE


class PharmCastMissing(RuntimeError):
    """Raised with instructions rather than a stack trace."""

    def __init__(self):
        super().__init__(
            "Cross-family reach needs PharmCast, which fingerprints the query "
            "the same way the published index was built.\n"
            "  pip install git+https://github.com/smuskal/pharmcast.git\n"
            "  https://pharmcast.ai   https://github.com/smuskal/pharmcast\n"
            "Every other part of toxpred works without it.")


def _pharmcast():
    try:
        import pharmcast
    except ImportError:
        raise PharmCastMissing()
    return pharmcast


def fetch_model(home: Path) -> Path:
    """The published PharmCast checkpoint, verified against its SHA256SUMS."""
    import hashlib

    import requests

    home = Path(home)
    dest = home / Path(MODEL_URL).name
    from .data import _headers
    sums = requests.get(SUMS_URL, timeout=120, headers=_headers()).text
    want = None
    published = []
    for line in sums.splitlines():
        parts = line.split()
        if len(parts) < 2 or line.lstrip().startswith("#"):
            continue
        name = parts[-1].lstrip("*./")
        published.append(name)
        if name.endswith(dest.name):
            want = parts[0]
    # An unlisted checkpoint is refused whether or not a copy is already
    # cached. Using a cached file that cannot be checked is the one outcome
    # worse than failing, because the run would look verified and not be.
    if want is None:
        raise RuntimeError(
            "%s is not listed in %s, so it cannot be verified.\n"
            "Published right now: %s"
            % (dest.name, SUMS_URL, ", ".join(published) or "nothing"))
    if not dest.exists():
        with requests.get(MODEL_URL, stream=True, timeout=1200,
                          headers=_headers()) as r:
            r.raise_for_status()
            with open(dest, "wb") as fh:
                for chunk in r.iter_content(1 << 20):
                    fh.write(chunk)
    if want:
        h = hashlib.sha256()
        with open(dest, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        if h.hexdigest() != want:
            dest.unlink(missing_ok=True)
            raise RuntimeError("PharmCast checkpoint failed its checksum")
    return dest


def targets_by_component(sites_tsv) -> dict:
    """component id -> the UniProt accessions it was solved against."""
    by_comp = collections.defaultdict(set)
    with open(sites_tsv) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ci, ai = head.index("component"), head.index("accession")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > max(ci, ai):
                by_comp[f[ci]].add(f[ai])
    return by_comp


def _data_file(name) -> Path:
    return Path(__file__).resolve().parent / "data_files" / name


def default_family_map() -> dict:
    """The level 2 family of each accession, as the paper counted them.

    This ships with the package, so families reached works with nothing
    supplied. It is the grouping of the released cross-family comparator at
    familyfoundationmodel.com, with a target carrying two families appearing
    twice and counting as two. Pass --family-map to group them your own way.
    """
    return read_family_map(_data_file("families.tsv"))


def default_antitargets() -> dict:
    """The 47 safety pharmacology antitargets, accession -> name.

    antitargets_hit counts how many of these the reached proteins include. It
    is the fourth of the four integers a contributor can release in place of a
    structure, and this is the list every number in the paper used.
    """
    out = {}
    with open(_data_file("antitargets.tsv")) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                out[parts[0].strip()] = parts[1].strip()
    return out


def read_family_map(path) -> dict:
    """Optional: a two column file, accession then family.

    Families are grouped however you choose to group them. Repeat an accession
    to put one protein in two families; a protein in two families counts as two
    rather than as a third family of its own.
    """
    out = collections.defaultdict(set)
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.replace(",", "\t").split("\t")
            if len(parts) >= 2:
                out[parts[0].strip()].add(parts[1].strip())
    return out


def reach(smiles, index_pfp, sites_tsv, model_path, family_map=None,
          cutoff: float = 0.5, block: int = 256, top_k: int = 500,
          antitargets=None):
    """-> one dict per query: how far its chemistry reaches.

    The operating point is the published one: the top_k most similar indexed
    ligands are taken first and then floored at cutoff, which is the order the
    paper's screen used. Raising top_k past 500 changes the count.

    Families reached uses the shipped grouping unless family_map is given, so
    the count is a count of families by default. Pass family_map=False to count
    only proteins. antitargets defaults to the shipped 47.
    """
    if family_map is None:
        family_map = default_family_map()
    elif family_map is False:
        family_map = None
    if antitargets is None:
        antitargets = default_antitargets()
    pc = _pharmcast()
    names, words = [], []
    for name, w in pc.read_pfp(index_pfp):
        names.append(name)
        words.append(w)
    tgts = targets_by_component(sites_tsv)
    comp_tgts = [sorted(tgts.get(n, ())) for n in names]

    model = pc.PharmCast.load(model_path)
    smiles = [smiles] if isinstance(smiles, str) else list(smiles)
    out = []
    for s in range(0, len(smiles), block):
        chunk = smiles[s:s + block]
        qwords = model.words_batch(chunk)
        T = pc.pharmtan_matrix(qwords, words)
        for i in range(len(chunk)):
            # The published order: most similar first, capped at top_k, then
            # floored at cutoff. Flooring first and capping after would take a
            # different 500 on a compound with more than 500 matches.
            row = T[i]
            order = sorted(range(len(row)), key=lambda j: -row[j])[:top_k]
            hit = [j for j in order if row[j] >= cutoff]
            reached = set()
            for j in hit:
                reached.update(comp_tgts[j])
            rec = dict(smiles=chunk[i],
                       targets_reached=len(reached),
                       matches=len(hit),
                       antitargets_reached=len(reached & set(antitargets)),
                       best_similarity=float(max(row)) if len(row) else 0.0)
            if family_map:
                fams = {f for a in reached for f in family_map.get(a, ())}
                rec["families_reached"] = len(fams)
                rec["families"] = sorted(fams)
            out.append(rec)
    return out
