"""Figure 6. Is the score good enough to use as a filter?

Panel A answers it in the only terms that matter to someone holding a virtual
library: out of 100 compounds, how many are toxic now, and how many are toxic in
the tenth the score flags. The gap is the whole value of the method.

Panel B is the catch against the cost. Discard the worst scoring compounds and the
line says what share of the toxic ones you removed. The diagonal is what you would
get by discarding at random, so distance above the diagonal is the gain.

Panel C sets the three lights and reports what lands in each.

    python code/make_figure6_useful.py
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# In the published package the measured result files this reads ship beside it
# in figures/data, and the figures are written to figures/out. Both can be
# pointed elsewhere, which is how the analysis tree runs the same script over a
# fresh run without editing it.
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.environ.get("TOXPRED_FIGURE_DATA") or os.path.join(HERE, "data")
OUTDIR = os.environ.get("TOXPRED_FIGURE_OUT") or os.path.join(HERE, "out")
os.makedirs(OUTDIR, exist_ok=True)
OUT = os.path.join(OUTDIR, "figure6_usefulness.png")

BEFORE, AFTER = "#8a8a8a", "#eb6834"
RED, AMBER, GREEN = "#c0392b", "#e8a33d", "#1baf7a"
INK, MUTED, GRID = "#1a1a19", "#52514e", "#d9d8d2"
CURVE_EP = ["CYP2D6", "CYP3A4", "Clinical toxicity", "Hepatotoxicity"]
BAND_EP = ["CYP2D6", "CYP3A4", "Clinical toxicity", "Hepatotoxicity"]


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
    d = json.load(open(os.path.join(RES, "phase12_usefulness.json")))
    S = sorted(d["summary"], key=lambda r: r["enrichment"])
    C = d["curves"]

    fig = plt.figure(figsize=(11.4, 10.2), dpi=200)
    fig.patch.set_facecolor("white")
    fig.text(0.055, 0.978, "Score a library, set the worst aside", fontsize=14,
             color=INK, ha="left", fontweight="bold")
    fig.text(0.055, 0.950,
             "34,197 compounds, 28 endpoints. Every number is measured on "
             "compounds the score was never fitted on.",
             fontsize=9.5, color=MUTED, ha="left", va="top")

    # ---- A: before and after ----
    a = fig.add_axes([0.26, 0.475, 0.70, 0.400])
    a.set_facecolor("white")
    y = np.arange(len(S))
    for i, r in enumerate(S):
        a.plot([r["percent_toxic_overall"], r["percent_toxic_in_flagged_tenth"]],
               [i, i], color=GRID, lw=1.8, zorder=1, solid_capstyle="round")
    a.scatter([r["percent_toxic_overall"] for r in S], y, s=30, color=BEFORE,
              zorder=3, edgecolor="white", linewidth=0.8,
              label="in the library as it stands")
    a.scatter([r["percent_toxic_in_flagged_tenth"] for r in S], y, s=36,
              color=AFTER, zorder=3, edgecolor="white", linewidth=0.8,
              label="in the tenth the score flags")
    a.set_yticks(y)
    a.set_yticklabels(["%s   (%.1fx)" % (r["endpoint"], r["enrichment"])
                       for r in S], fontsize=8)
    a.set_xlabel("out of every 100 compounds, how many are toxic", fontsize=9.5)
    a.set_title("A   Score the library, set aside the worst tenth", fontsize=10,
                color=INK, loc="left", pad=26)
    a.legend(fontsize=8.5, frameon=False, loc="upper center",
             bbox_to_anchor=(0.5, 1.085), ncol=2)
    a.set_xlim(0, 105)
    a.set_xticks([0, 20, 40, 60, 80, 100])
    a.set_ylim(-0.8, len(S) - 0.2)
    style(a, grid="x")

    # ---- B: catch against cost ----
    b = fig.add_axes([0.075, 0.215, 0.38, 0.195])
    b.set_facecolor("white")
    cmap = plt.get_cmap("tab10")
    b.plot([0, 100], [0, 100], color=MUTED, lw=1.3, ls=(0, (4, 3)), zorder=2)
    b.annotate("discarding at random", xy=(58, 52), fontsize=8, color=MUTED,
               rotation=32)
    for i, ep in enumerate(CURVE_EP):
        if ep not in C:
            continue
        x = 100 * np.asarray(C[ep]["fraction_flagged"])
        yy = 100 * np.asarray(C[ep]["toxics_caught"])
        b.plot(x, yy, "-", lw=2.2, color=cmap(i), zorder=3, label=ep)
    b.set_xlabel("percent of the library discarded", fontsize=9)
    b.set_ylabel("percent of the toxic\ncompounds removed", fontsize=9)
    b.set_title("B   How much you catch for what you give up", fontsize=10, color=INK, loc="left",
                pad=6)
    b.legend(fontsize=8, frameon=False, loc="lower right")
    b.set_xlim(0, 100)
    b.set_ylim(0, 102)
    style(b)

    # ---- C: the three lights ----
    c = fig.add_axes([0.575, 0.215, 0.385, 0.195])
    c.set_facecolor("white")
    eps = [e for e in BAND_EP if any(r["endpoint"] == e for r in S)]
    w = 0.26
    for j, ep in enumerate(eps):
        r = [x for x in S if x["endpoint"] == ep][0]
        for k, (band, colour) in enumerate(zip(r["bands"], (RED, AMBER, GREEN))):
            c.bar(j + (k - 1) * w, band["percent_toxic"], width=w * 0.92,
                  color=colour, zorder=3,
                  label=band["band"].split(",")[0] if j == 0 else None)
        c.plot([j - 1.6 * w, j + 1.6 * w],
               [r["percent_toxic_overall"]] * 2, color=INK, lw=1.4, zorder=4)
    c.annotate("black line: the library\nbefore any filtering", xy=(0.02, 0.86),
               xycoords="axes fraction", fontsize=8, color=MUTED, va="top")
    c.set_xticks(range(len(eps)))
    c.set_xticklabels(eps, fontsize=8.5)
    c.set_ylabel("percent toxic in the band", fontsize=9)
    c.set_title("C   A red, amber and green split", fontsize=10, color=INK, loc="left",
                pad=6)
    c.legend(fontsize=8, frameon=False, loc="upper right", ncol=3,
             columnspacing=0.9, handlelength=1.1)
    c.set_ylim(0, 108)
    c.set_yticks([0, 25, 50, 75, 100])
    c.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])
    style(c)

    # Read the examples out of the results rather than typing them in. A
    # hardcoded "45 in 100" in this caption survived a rerun that moved it to
    # 43, which is how a figure starts disagreeing with its own data.
    by_ep = {r["endpoint"]: r for r in d["summary"]}

    def pair(name):
        r = by_ep[name]
        return (round(r["percent_toxic_overall"]),
                round(r["percent_toxic_in_flagged_tenth"]))

    d6, d6f = pair("CYP2D6")
    ct, ctf = pair("Clinical toxicity")
    hep = round(by_ep["Hepatotoxicity"]["percent_toxic_overall"])
    c30 = round(by_ep["Cardiotoxicity-30"]["percent_toxic_overall"])
    fig.text(0.055, 0.145,
             "How to read this. The multiplier beside each endpoint in panel A is "
             "how much richer the flagged tenth is than the library it came from.\n"
             "The method earns its place where that number is large AND the "
             "starting rate is low: CYP2D6 goes from %d compounds in 100 to %d in "
             "100,\nand clinical toxicity from %d in 100 to %d in 100. It earns "
             "little where most of the library is already toxic, which is why\n"
             "hepatotoxicity at %d percent and cardiotoxicity at %d percent barely "
             "move. A filter cannot enrich what is already everywhere."
             % (d6, d6f, ct, ctf, hep, c30),
             fontsize=9, color=INK, ha="left", va="top")

    fig.savefig(OUT, facecolor="white")
    print("wrote %s" % OUT)


if __name__ == "__main__":
    main()
