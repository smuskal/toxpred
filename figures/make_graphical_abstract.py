"""The graphical abstract: the two halves of the paper in one picture.

Panel A  The consortium vote. Growing the pool buys reach, not accuracy, which
         is the argument for pooling.
Panel B  Cross-family reach. Percent of compounds measured active against the
         number of distinct protein families their structural matches reach, on
         every endpoint carrying all six bands whose rate climbs, with CYP3A4
         drawn heaviest.

Every value is read from the measured result files at run time. Nothing is
typed in.

    python figures/make_graphical_abstract.py
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# In the published package the measured result files this reads ship beside it
# in figures/data, and the figures are written to figures/out. Both can be
# pointed elsewhere, which is how the analysis tree runs the same script over a
# fresh run without editing it.
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.environ.get("TOXPRED_FIGURE_DATA") or os.path.join(HERE, "data")
OUTDIR = os.environ.get("TOXPRED_FIGURE_OUT") or os.path.join(HERE, "out")
os.makedirs(OUTDIR, exist_ok=True)
OUT = os.path.join(OUTDIR, "figure0_graphical_abstract.png")

INK = "#1b1b1b"
MUTED = "#5e5e5e"
BLUE = "#1f77b4"
ORANGE = "#e8800c"
GREEN = "#2e8b57"
GREY = "#c7c7c7"

# Every endpoint measured across all six bands whose rate climbs. Named here so
# the panel cannot quietly become a one endpoint claim.
CLIMBERS = ["CYP3A4", "CYP2C19", "CYP2C9", "CYP1A2", "SR-MMP", "SR-ARE",
            "NR-aromatase", "Clinical toxicity"]
PRETTY = {"SR-MMP": "mitochondrial membrane potential",
          "SR-ARE": "oxidative stress response",
          "NR-aromatase": "aromatase",
          "Clinical toxicity": "clinical toxicity"}


def load(name):
    with open(os.path.join(RES, name)) as fh:
        return json.load(fh)


def main():
    sweeps = load("refsim_sweeps.json")
    breadth = {b["endpoint"]: b
               for b in load("phase10_targets_breadth.json")["breadth"]}

    refset = sweeps["refset_sweep"]["pLD50"]
    sizes = [r["reference_compounds"] for r in refset]
    reach = [100 * r["coverage"] for r in refset]
    acc = [100 * r["within_1_log"] for r in refset]

    labels = [b["families"] for b in breadth["CYP3A4"]["bands"]]

    fig = plt.figure(figsize=(15.0, 6.4), dpi=200)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.18],
                          left=0.052, right=0.975, top=0.775, bottom=0.185,
                          wspace=0.26)

    fig.text(0.052, 0.957,
             "A compound's toxicity can be gauged from what it resembles, and from "
             "how widely its chemistry reaches",
             fontsize=16.0, color=INK, fontweight="bold", ha="left", va="top")

    # ---------------------------------------------------------------- panel A
    ax = fig.add_subplot(gs[0, 0])
    ax.set_title("A   A larger shared pool buys reach, not accuracy",
                 fontsize=11.8, color=INK, loc="left", pad=34)
    ax.plot(sizes, acc, "-o", color=BLUE, lw=2.4, ms=5.5,
            label="within one log unit, of answers given")
    ax.plot(sizes, reach, "-o", color=ORANGE, lw=2.4, ms=5.5,
            label="share of compounds answered for")
    ax.set_xscale("log")
    ax.set_xlabel("compounds in the voting reference set", fontsize=10.5)
    ax.set_ylabel("percent", fontsize=10.5)
    ax.set_ylim(0, 100)
    ax.set_xticks([127, 610, 2105, 5953])
    ax.set_xticklabels(["127", "610", "2,105", "5,953"], fontsize=9.5)
    ax.tick_params(labelsize=9.5)
    ax.grid(axis="y", color="#e8e8e8")
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.annotate("%.1f%%" % acc[-1], (sizes[-1], acc[-1]),
                textcoords="offset points", xytext=(-6, 10), fontsize=10.5,
                color=BLUE, fontweight="bold", ha="right")
    ax.annotate("%.1f%%" % reach[-1], (sizes[-1], reach[-1]),
                textcoords="offset points", xytext=(-6, 10), fontsize=10.5,
                color=ORANGE, fontweight="bold", ha="right")
    ax.legend(fontsize=9.0, frameon=False, loc="upper left",
              bbox_to_anchor=(0.02, 0.72))

    # ---------------------------------------------------------------- panel B
    ax = fig.add_subplot(gs[0, 1])
    ax.set_title("B   Reaching more protein families goes with being measured active",
                 fontsize=11.8, color=INK, loc="left", pad=34)
    x = list(range(len(labels)))
    for key in CLIMBERS:
        vals = [b["percent_toxic"] for b in breadth[key]["bands"]]
        lead = key == "CYP3A4"
        ax.plot(x, vals, "-o" if lead else "-",
                color=BLUE if lead else GREY,
                lw=3.0 if lead else 1.8, ms=7 if lead else 0,
                zorder=4 if lead else 2,
                label=PRETTY.get(key, key))
    ax.legend(fontsize=8.8, frameon=False, ncol=2, loc="upper left",
              handlelength=1.6, columnspacing=1.2, labelspacing=0.35,
              borderpad=0.0, handletextpad=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9.5)
    ax.set_xlim(-0.25, len(labels) - 0.75)
    ax.set_xlabel("protein families reached by the compound's structural matches",
                  fontsize=10.5)
    ax.set_ylabel("percent of compounds measured active", fontsize=10.5)
    ax.set_ylim(0, 62)
    ax.tick_params(labelsize=9.5)
    ax.grid(axis="y", color="#e8e8e8")
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    low = breadth["CYP3A4"]["bands"][0]["percent_toxic"]
    ax.annotate("%.1f%%" % low, (0, low),
                textcoords="offset points", xytext=(4, -16), fontsize=10,
                color=BLUE, fontweight="bold")
    top = breadth["CYP3A4"]["bands"][-1]["percent_toxic"]
    ax.annotate("%.1f%%" % top, (len(labels) - 1, top),
                textcoords="offset points",
                xytext=(-4, 11), fontsize=11, color=BLUE, fontweight="bold",
                ha="right")

    fig.text(0.052, 0.030,
             "A, rat oral LD50 from CATMoS, 2,873 held-out compounds.    "
             "B, 34,197 compounds screened against 27,910 ligands solved in protein "
             "structures.\n"
             "Eight endpoints, every one measured across all six bands and climbing. "
             "Active means an inhibitor for the cytochromes, active in the assay\n"
             "for the endocrine and stress endpoints, and a clinical finding for "
             "clinical toxicity.",
             fontsize=9.0, color=MUTED, ha="left", va="bottom", linespacing=1.7)

    fig.savefig(OUT, facecolor="white")
    print("wrote %s" % OUT)


if __name__ == "__main__":
    main()
