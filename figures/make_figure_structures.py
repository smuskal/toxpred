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
import os
import textwrap

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

INK, MUTED = "#1a1a19", "#52514e"
COLS, ROWS = 3, 4
PANEL = (520, 380)


def rows():
    with open(PROV) as fh:
        got = list(csv.DictReader(fh))
    if len(got) != COLS * ROWS:
        raise SystemExit("%s holds %d compounds and the figure is %d panels"
                         % (PROV, len(got), COLS * ROWS))
    return got


def draw(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise SystemExit("this SMILES does not parse: %s" % smiles)
    return Draw.MolToImage(m, size=PANEL, kekulize=True)


def main():
    got = rows()
    fig = plt.figure(figsize=(13.5, 16.5), dpi=150)
    # Two axes per panel: the structure above, the text below. Nothing is drawn
    # outside its own axis, which is what keeps the blocks from overlapping.
    gs = fig.add_gridspec(ROWS * 2, COLS, height_ratios=[3, 2] * ROWS,
                          hspace=0.08, wspace=0.06, top=0.965, bottom=0.01,
                          left=0.01, right=0.99)
    for i, r in enumerate(got):
        row, col = divmod(i, COLS)
        ax = fig.add_subplot(gs[row * 2, col])
        ax.imshow(draw(r["smiles"]))
        ax.axis("off")

        tx = fig.add_subplot(gs[row * 2 + 1, col])
        tx.axis("off")
        toxic = [e for e in (r["endpoints_toxic"] or "").split("|") if e]
        name = r["name_chembl37"]
        head = "%s\n%s" % (name, r["inchikey"])
        body = ("%s families, %s proteins reached\ntoxic on %s of %s measured\n%s"
                % (r["families_reached"], r["proteins_reached"], r["n_toxic"],
                   r["n_measured"], textwrap.fill(", ".join(toxic), 44)))
        tx.text(0, 1, head, transform=tx.transAxes, va="top", ha="left",
                fontsize=10.5, color=INK, fontweight="bold", linespacing=1.5)
        tx.text(0, 0.62, body, transform=tx.transAxes, va="top", ha="left",
                fontsize=9.5, color=MUTED, linespacing=1.6)

    fig.suptitle("Twelve compounds in the top band of families reached",
                 fontsize=13, color=INK, y=0.995)
    fig.savefig(OUT, facecolor="white")
    print("wrote %s" % OUT)


if __name__ == "__main__":
    main()
