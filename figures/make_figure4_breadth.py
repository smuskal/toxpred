"""Figure 4. Reaching more protein families goes with being more often toxic.

Panel A is the whole claim in one picture: group every compound by how many
distinct protein families its reverse screen reaches in the crystallographic
record, then read off what percentage of each group is toxic. No model, no
baseline, no statistic. Counts and percentages.

Panel B names the proteins. For one endpoint, the proteins whose co-crystallised
ligands most resemble the toxic compounds, with the percentage of toxic and of
non toxic compounds that reach each one.

    python code/make_figure4_breadth.py
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def index_facts():
    """The release every figure here was drawn against."""
    return json.load(open(os.path.join(RES, "index_release.json")))


def target_names():
    """Accession -> short protein name, for the human proteins panel C labels.

    The published target table restricts this to human, so an accession absent
    from it is left unlabeled and drops out of the panel, which is what the
    panel did when it read that table directly.
    """
    out = {}
    with open(os.path.join(RES, "target_names.tsv")) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) >= 2:
                out[f[0].strip()] = f[1].strip()
    return out

# In the published package the measured result files this reads ship beside it
# in figures/data, and the figures are written to figures/out. Both can be
# pointed elsewhere, which is how the analysis tree runs the same script over a
# fresh run without editing it.
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.environ.get("TOXPRED_FIGURE_DATA") or os.path.join(HERE, "data")
OUTDIR = os.environ.get("TOXPRED_FIGURE_OUT") or os.path.join(HERE, "out")
os.makedirs(OUTDIR, exist_ok=True)
OUT = os.path.join(OUTDIR, "figure4_breadth_and_targets.png")

TOXIC, CLEAN = "#eb6834", "#2a78d6"
INK, MUTED, GRID = "#1a1a19", "#52514e", "#d9d8d2"
SHOW = ["CYP3A4", "CYP2C9", "CYP2C19", "CYP1A2", "Cardiotoxicity-5",
        "Clinical toxicity"]
DETAIL = "CYP3A4"


def style(ax, grid="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.grid(axis=grid, color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=8.5, colors=MUTED, length=0)
    ax.xaxis.label.set_color(MUTED)
    ax.yaxis.label.set_color(MUTED)


def main():
    d = json.load(open(os.path.join(RES, "phase10_targets_breadth.json")))
    by_ep = {r["endpoint"]: r for r in d["breadth"]}
    prot = {r["endpoint"]: r for r in d["proteins"]}

    fig = plt.figure(figsize=(11.0, 9.6), dpi=200)
    fig.patch.set_facecolor("white")
    fig.text(0.06, 0.968, "Reaching more protein families goes with being more "
             "often toxic", fontsize=13, color=INK, ha="left",
             fontweight="bold")
    fig.text(0.06, 0.944,
             "Every compound is screened against %s ligands solved in protein "
             "structures.\nFamilies reached is how many distinct protein families "
             "those matches belong to, at a required similarity of 0.50."
             % format(index_facts()["ligands"], ","),
             fontsize=9, color=MUTED, ha="left", va="top")

    a = fig.add_axes([0.075, 0.585, 0.88, 0.315])
    a.set_facecolor("white")
    shown = [e for e in SHOW if e in by_ep]
    cmap = plt.get_cmap("tab10")
    ends = []
    for i, ep in enumerate(shown):
        r = by_ep[ep]
        x = [b["low"] for b in r["bands"]]
        y = [b["percent_toxic"] for b in r["bands"]]
        a.plot(x, y, "-o", lw=2.2, ms=6, color=cmap(i), zorder=3,
               markeredgecolor="white", markeredgewidth=1.1)
        ends.append([y[-1], x[-1], ep, cmap(i)])
    # push end labels apart so none can overlap, keeping their order
    ends.sort(key=lambda e: e[0])
    gap = 3.4
    for j in range(1, len(ends)):
        if ends[j][0] - ends[j - 1][0] < gap:
            ends[j][0] = ends[j - 1][0] + gap
    for lab_y, x_end, ep, colour in ends:
        a.annotate(" %s" % ep, xy=(x_end, lab_y), fontsize=9, color=colour,
                   va="center", fontweight="bold")
    a.set_xticks([b["low"] for b in by_ep[shown[0]]["bands"]])
    a.set_xticklabels([b["families"] for b in by_ep[shown[0]]["bands"]],
                      fontsize=8.5)
    a.set_xlabel("protein families reached by the compound's structural matches",
                 fontsize=9.5)
    a.set_ylabel("percent of compounds that are toxic", fontsize=9.5)
    a.set_title("A", fontsize=10, color=INK, loc="left", pad=6)
    a.set_xlim(-0.8, 23.5)
    a.set_yticks([0, 20, 40, 60, 80])
    a.set_yticklabels(["0%", "20%", "40%", "60%", "80%"])
    style(a)

    b = fig.add_axes([0.30, 0.135, 0.655, 0.375])
    b.set_facecolor("white")
    r = prot[DETAIL]
    human = target_names()
    for t in r["top"]:
        t["protein"] = human.get(t["accession"], "")
    top = [t for t in r["top"] if t["protein"]][:10]
    top = list(reversed(top))
    y = range(len(top))
    for i, t in enumerate(top):
        b.plot([t["percent_of_non_toxic"], t["percent_of_toxic"]], [i, i],
               color=GRID, lw=1.8, zorder=1, solid_capstyle="round")
    b.scatter([t["percent_of_non_toxic"] for t in top], y, s=34, color=CLEAN,
              marker="o", zorder=3, edgecolor="white", linewidth=0.8,
              label="of compounds that are NOT %s inhibitors" % DETAIL)
    b.scatter([t["percent_of_toxic"] for t in top], y, s=38, color=TOXIC,
              marker="o", zorder=3, edgecolor="white", linewidth=0.8,
              label="of compounds that ARE %s inhibitors" % DETAIL)
    b.set_yticks(list(y))
    b.set_yticklabels([t["protein"][:42] for t in top], fontsize=8.5)
    b.set_xlabel("percent of compounds whose structural matches reach that protein",
                 fontsize=9.5)
    b.set_title("B   The proteins most associated with %s inhibition, named"
                % DETAIL, fontsize=10, color=INK, loc="left", pad=6)
    b.legend(fontsize=8.5, frameon=False, loc="upper center",
             bbox_to_anchor=(0.5, -0.145), ncol=2)
    b.set_xlim(0, 100)
    b.set_xticks([0, 20, 40, 60, 80, 100])
    b.set_xticklabels(["0%", "20%", "40%", "60%", "80%", "100%"])
    b.set_ylim(-0.8, len(top) - 0.2)
    style(b, grid="x")

    fig.text(0.06, 0.028,
             "Panel B reads: of the compounds that inhibit CYP3A4, this percentage "
             "have a structural match sitting in that protein. These are the sites "
             "whose\ncrystallised ligands most resemble the inhibitors. Whether "
             "these proteins carry the toxicity is a separate question.",
             fontsize=8.5, color=MUTED, ha="left", va="bottom")

    fig.savefig(OUT, facecolor="white")
    print("wrote %s" % OUT)


if __name__ == "__main__":
    main()
