"""Twelve measured compounds from the top band of the reach figure.

Every field is read from the record the figure was built from,
data/figure_structures_provenance.csv, and nothing is typed in here:

    name                 ChEMBL 37 preferred name, joined by standard InChIKey
    smiles               the structure the panel draws and the InChIKey it was
                         selected by
    families, proteins   the reverse screen at a required similarity of 0.50
                         against the published index
    endpoints            the nine binding-mechanism TOXRIC tables, as the
                         endpoints measured and the endpoints toxic

THE DENOMINATOR IS THE POINT. A compound toxic on three of the three endpoints
it was tested on is not making the same claim as one toxic on three of nine, so
every panel prints both. Without it a broadly toxic compound cannot be told
apart from a rarely tested one.

Each panel owns a structure axis and a text axis, and text is drawn only inside
its own, so a long endpoint list cannot run into the panel below.

    python figures/make_figure_structures.py
"""
import csv
import json
import os
import textwrap

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rdkit import Chem, RDLogger
from rdkit.Chem import Draw

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.environ.get("TOXPRED_FIGURE_DATA") or os.path.join(HERE, "data")
OUTDIR = os.environ.get("TOXPRED_FIGURE_OUT") or os.path.join(HERE, "out")
os.makedirs(OUTDIR, exist_ok=True)
PROV = os.path.join(RES, "figure_structures_provenance.csv")
OUT = os.path.join(OUTDIR, "figure_structures.png")

INK = "#171a1e"
REACH = "#0b5f77"
TOXIC = "#a33227"
MUTED = "#6b7280"
NR, NC = 3, 4


def rows():
    with open(PROV) as fh:
        got = list(csv.DictReader(fh))
    if len(got) != NR * NC:
        raise SystemExit("%s holds %d compounds and the figure is %d panels"
                         % (PROV, len(got), NR * NC))
    return got


def mol(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise SystemExit("this SMILES does not parse: %s" % smiles)
    return m


def main():
    got = rows()
    fig = plt.figure(figsize=(15.0, 13.4), dpi=200)
    fig.patch.set_facecolor("white")
    # Two rows of axes per compound row: the drawing, then its own text box.
    # Nothing is drawn outside its own axis, which is what keeps the blocks
    # from overlapping.
    gs = fig.add_gridspec(NR * 2, NC, height_ratios=[1.0, 0.80] * NR,
                          top=0.925, bottom=0.062, left=0.012, right=0.988,
                          hspace=0.10, wspace=0.04)

    for i, r in enumerate(got):
        row, col = divmod(i, NC)
        ax = fig.add_subplot(gs[row * 2, col])
        ax.imshow(np.asarray(Draw.MolToImage(mol(r["smiles"]), size=(760, 480))))
        ax.axis("off")
        ax.set_title(r["name_chembl37"], fontsize=15.5, fontweight="bold",
                     color=INK, pad=4)

        toxic = [e for e in (r["endpoints_toxic"] or "").split("|") if e]
        tx = fig.add_subplot(gs[row * 2 + 1, col])
        tx.axis("off")
        tx.set_xlim(0, 1)
        tx.set_ylim(0, 1)
        tx.text(0.5, 1.00, "reaches %s families, %s proteins"
                % (r["families_reached"], r["proteins_reached"]),
                ha="center", va="top", fontsize=11.6, color=REACH,
                fontweight="semibold")
        tx.text(0.5, 0.78, "active on %s of %s endpoints"
                % (r["n_toxic"], r["n_measured"]),
                ha="center", va="top", fontsize=11.6, color=TOXIC,
                fontweight="semibold")
        tx.text(0.5, 0.56, textwrap.fill(", ".join(toxic), 40),
                ha="center", va="top", fontsize=9.8, color=TOXIC,
                linespacing=1.3)
        # fixed offset, not bottom-anchored: a bottom anchor let this line
        # spill out of its own axis and into the next row's compound name.
        tx.text(0.5, 0.17, r["inchikey"], ha="center", va="top", fontsize=7.4,
                color=MUTED, family="DejaVu Sans Mono")

    fig.suptitle("Compounds reaching many protein families are measured active "
                 "on many of the nine assays", fontsize=20, fontweight="bold",
                 color=INK, y=0.972)

    facts = json.load(open(os.path.join(RES, "index_release.json")))
    caption = (
        "Twelve of the 34,197 compounds in the TOXRIC tables, all of them in the "
        "top band of Figure 3. Families and proteins reached are counted by the "
        "reverse screen at a required similarity of 0.50 against the published "
        "index of %s ligands solved in protein structures, release "
        "%s. Nine endpoints are counted: CYP3A4, CYP2C9, CYP2C19 and CYP1A2 "
        "inhibition, NR-AR and NR-ER binding, Ames mutagenicity, hepatotoxicity, "
        "and cardiotoxicity at 1 micromolar. The denominator is how many of those "
        "nine carry a measurement for that compound. Names are "
        "ChEMBL 37 preferred names and each panel carries the InChIKey it was "
        "selected by." % (format(facts["ligands"], ","), facts["release"]))
    fig.text(0.5, 0.012, textwrap.fill(caption, 148), ha="center", va="bottom",
             fontsize=9.6, color="#444444", linespacing=1.45)

    fig.savefig(OUT, facecolor="white")
    plt.close(fig)
    print("wrote %s" % OUT)


if __name__ == "__main__":
    main()
