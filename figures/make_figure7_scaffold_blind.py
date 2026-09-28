"""Figure 7. The trend survives forbidding the screen its own analogs.

Three screens per compound: every indexed ligand available, every ligand sharing
the query's scaffold removed, and every ligand at or above 0.70 similarity
removed. The three lines lie on top of each other, which is the result.
"""
import json, os
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
OUT = os.path.join(OUTDIR, "figure7_scaffold_blind.png")
INK, MUTED, GRID = "#1a1a19", "#52514e", "#d9d8d2"
SHOW = ["CYP3A4", "CYP2C9", "CYP2C19", "Clinical toxicity"]
MODES = [("as before", "-", 2.6, 1.0),
         ("scaffold blind", "--", 2.0, 0.95),
         ("analogue blind", ":", 2.0, 0.95)]
# The mode strings above are keys into the result file and must stay as they
# are written there. These are what the legend prints.
LABEL = {"analogue blind": "analog blind"}


def main():
    d = json.load(open(os.path.join(RES, "phase13_scaffold_blind.json")))
    by = {r["endpoint"]: r for r in d["endpoints"]}
    fig, axes = plt.subplots(1, 4, figsize=(13.2, 4.0), dpi=200)
    fig.patch.set_facecolor("white")
    cmap = plt.get_cmap("tab10")
    for k, ep in enumerate(SHOW):
        ax = axes[k]
        ax.set_facecolor("white")
        r = by[ep]
        for mi, (mode, ls, lw, al) in enumerate(MODES):
            b = r["modes"][mode]
            ax.plot([x["low"] for x in b], [x["percent_toxic"] for x in b],
                    ls, lw=lw, alpha=al, color=cmap(k), zorder=3 - mi,
                    marker="o" if mi == 0 else None, ms=5,
                    markeredgecolor="white", markeredgewidth=1.0,
                    label=LABEL.get(mode, mode) if k == 0 else None)
        ax.set_title("%s\n%.1f%% toxic overall" % (ep, r["percent_toxic_overall"]),
                     fontsize=9.5, color=INK, loc="left", pad=6)
        ax.set_xticks([0, 4, 8, 12, 16])
        ax.set_xticklabels(["0", "4-7", "8-11", "12-15", "16+"], fontsize=8)
        ax.set_ylim(0, 70)
        ax.set_yticks([0, 20, 40, 60])
        ax.set_yticklabels(["0%", "20%", "40%", "60%"], fontsize=8.5)
        if k == 0:
            ax.set_ylabel("percent of compounds\nthat are toxic", fontsize=9)
            ax.legend(fontsize=8.5, frameon=False, loc="upper left")
        ax.set_xlabel("families reached", fontsize=9)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(GRID)
        ax.grid(axis="y", color=GRID, lw=0.8, zorder=0)
        ax.set_axisbelow(True)
        ax.tick_params(labelsize=8.5, colors=MUTED, length=0)
        ax.xaxis.label.set_color(MUTED); ax.yaxis.label.set_color(MUTED)
    fig.suptitle("The trend is not the compound finding its own analogs",
                 fontsize=13, fontweight="bold", color=INK, x=0.012, ha="left",
                 y=0.985)
    # Read these out of the result file. They were typed in once and then went
    # stale by two to three points when the screen was rerun.
    eb = d.get("effect_of_blinding", {})
    sb = eb.get("scaffold blind", {})
    ab = eb.get("analogue blind", {})
    fig.text(0.012, 0.915,
             "%.1f percent of compounds share a scaffold with an indexed ligand "
             "and %.1f percent have a match at or above 0.70 similarity. Removing "
             "them changes\nfamilies reached for %.1f and %.1f percent of compounds "
             "respectively, by %.3f and %.3f families on average, and moves none of "
             "these curves."
             % (d.get("percent_sharing_a_scaffold", float("nan")),
                d.get("percent_with_a_match_at_or_above_0.70", float("nan")),
                sb.get("percent_of_compounds_changed", float("nan")),
                ab.get("percent_of_compounds_changed", float("nan")),
                sb.get("mean_change_in_families", float("nan")),
                ab.get("mean_change_in_families", float("nan"))),
             fontsize=9, color=MUTED, ha="left", va="top")
    fig.subplots_adjust(left=0.065, right=0.99, top=0.70, bottom=0.14, wspace=0.28)
    fig.savefig(OUT, facecolor="white")
    print("wrote %s" % OUT)


if __name__ == "__main__":
    main()
