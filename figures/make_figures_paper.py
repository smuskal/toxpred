"""Paper figures for the consortium work.

Figure 2 reproduces the 2003 paper's three sweeps and its error distribution on
public data, so the method is shown to behave as published before anything new is
claimed. Figure 3 is the new part: the same machinery over 28 further toxicity
endpoints, and the traffic light that comes out of it.

Accuracy and coverage are both proportions, so they share one axis and no second
scale is introduced.

    python code/make_figures_paper.py
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# In the published package the measured result files this reads ship beside it
# in figures/data, and the figures are written to figures/out. Both can be
# pointed elsewhere, which is how the analysis tree runs the same script over a
# fresh run without editing it.
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.environ.get("TOXPRED_FIGURE_DATA") or os.path.join(HERE, "data")
OUTDIR = os.environ.get("TOXPRED_FIGURE_OUT") or os.path.join(HERE, "out")
os.makedirs(OUTDIR, exist_ok=True)

PERF, COV = "#2a78d6", "#eb6834"
INK, MUTED, GRID = "#1a1a19", "#52514e", "#d9d8d2"


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


def two_series(ax, x, perf, cov, xlabel, title, perf_label, logx=False):
    perf = [p * 100 if p is not None and p == p else np.nan for p in perf]
    cov = [c * 100 for c in cov]
    ax.plot(x, perf, "-o", color=PERF, lw=2, ms=5.5, label=perf_label,
            markeredgecolor="white", markeredgewidth=1.1, zorder=3)
    ax.plot(x, cov, "-s", color=COV, lw=2, ms=5, label="compounds it will answer for",
            markeredgecolor="white", markeredgewidth=1.1, zorder=3)
    if logx:
        ax.set_xscale("log")
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_title(title, fontsize=9.5, color=INK, loc="left", pad=6)
    ax.set_ylim(0, 103)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.set_yticklabels(["0%", "20%", "40%", "60%", "80%", "100%"])
    style(ax)


def fig_2003(res):
    fig = plt.figure(figsize=(10.4, 7.6), dpi=200)
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.22, left=0.08,
                          right=0.975, top=0.845, bottom=0.085)
    fig.text(0.08, 0.972,
             "The consortium method on public data: what it answers, and how often",
             fontsize=12, color=INK, ha="left", fontweight="bold")
    fig.text(0.08, 0.928,
             "Rat oral LD50, CATMoS. Reference set 8,599 compounds, held-out set 2,873. "
             "Morgan radius 2, unweighted votes.\nDemanding a closer match trades "
             "coverage for accuracy. Growing the reference set improves both.",
             fontsize=9, color=MUTED, ha="left")

    a = fig.add_subplot(gs[0, 0])
    s = res["cutoff_sweep"]["pLD50"]
    two_series(a, [r["cutoff"] for r in s], [r["within_1_log"] for r in s],
               [r["coverage"] for r in s],
               "required similarity to a reference compound (Tanimoto)",
               "A   Demanding a closer match",
               "answers within one log unit of the measurement")
    a.legend(fontsize=8.5, frameon=False, loc="lower left")

    b = fig.add_subplot(gs[0, 1])
    s = res["members_sweep"]["pLD50"]
    x = [r["min_members"] for r in s]
    two_series(b, x, [r["within_1_log"] if r["within_1_log"] else np.nan for r in s],
               [r["coverage"] for r in s], "minimum consortium size",
               "B   Requiring more compounds to vote", "answers within one log unit of the measurement",
               logx=True)
    b.set_xticks([1, 2, 5, 10, 20, 50])
    b.set_xticklabels(["1", "2", "5", "10", "20", "50"])

    c = fig.add_subplot(gs[1, 0])
    s = res["refset_sweep"]["pLD50"]
    two_series(c, [r["reference_compounds"] for r in s],
               [r["within_1_log"] for r in s], [r["coverage"] for r in s],
               "reference compounds carrying an exact value",
               "C   Growing the shared reference set",
               "answers within one log unit of the measurement", logx=True)
    c.set_xticks([200, 500, 1000, 2000, 5000, 8600])
    c.set_xticklabels(["200", "500", "1,000", "2,000", "5,000", "8,600"])
    c.minorticks_off()

    d = fig.add_subplot(gs[1, 1])
    err = np.asarray(res["error_distribution"]["errors"])
    d.hist(err, bins=np.arange(-3, 3.01, 0.25), color=PERF, edgecolor="white",
           linewidth=0.8, zorder=3)
    d.axvline(0, color=MUTED, lw=1.2, zorder=2)
    within = float((np.abs(err) <= 1).mean())
    d.set_title("D   How wrong it is when it answers", fontsize=9.5, color=INK,
                loc="left", pad=6)
    d.set_xlabel("predicted minus observed pLD50, log$_{10}$ mg/kg", fontsize=9)
    d.annotate("%.1f%% within one log unit\n%s compounds answered at a required\n"
               "similarity of %.2f"
               % (100 * within, format(res["error_distribution"]["n"], ","),
                  res["error_distribution"]["cutoff"]),
               xy=(0.03, 0.78), xycoords="axes fraction", fontsize=8.5,
               color=MUTED)
    style(d)
    out = os.path.join(OUTDIR, "figure2_consortium_sweeps.png")
    fig.savefig(out, facecolor="white")
    print("wrote %s" % out)


def fig_endpoints(res):
    eps = res["endpoints"]
    rows = []
    for e in eps:
        for s in e["sweep"]:
            # Every endpoint at every cutoff. Dropping the points below 5%
            # coverage took endpoints out of the median at the tightest
            # matches, so the blue and orange lines were a median over 25
            # endpoints where the caption and the text say 28, and read 96.5%
            # and 9.0% at a cutoff of 0.90 against the 96.6% and 8.6%
            # reported. The operating point below is unaffected; it takes the
            # tightest cutoff still answering a fifth of compounds.
            if s["accuracy"] is None:
                continue
            rows.append(dict(category=e["category"], endpoint=e["endpoint"],
                             cutoff=s["cutoff"], accuracy=s["accuracy"],
                             coverage=s["coverage"],
                             caught=s.get("toxic_caught"),
                             passed=s.get("harmless_passed"),
                             balanced=s.get("balanced"),
                             mcc=s.get("mcc"),
                             n_toxic=s.get("n_toxic")))
    d = pd.DataFrame(rows)
    # operating point: the most demanding cutoff still answering a fifth of compounds
    op = (d[d.coverage >= 0.20].sort_values("cutoff")
          .groupby("endpoint").tail(1).reset_index(drop=True))
    op = op.dropna(subset=["mcc"]).sort_values("mcc").reset_index(drop=True)

    fig = plt.figure(figsize=(10.4, 9.8), dpi=200)
    fig.patch.set_facecolor("white")
    fig.text(0.06, 0.963,
             "The same machinery across 28 further toxicity endpoints",
             fontsize=12.5, color=INK, ha="left", fontweight="bold")
    fig.text(0.06, 0.936,
             "TOXRIC, 149,224 measurements over 34,162 compounds.",
             fontsize=9, color=MUTED, ha="left")

    a = fig.add_axes([0.085, 0.635, 0.87, 0.255])
    a.set_facecolor("white")
    cuts = sorted(d.cutoff.unique())
    acc = [100 * d.accuracy[d.cutoff == c].median() for c in cuts]
    cov = [100 * d.coverage[d.cutoff == c].median() for c in cuts]
    for ep, g in d.groupby("endpoint"):
        g = g.sort_values("cutoff")
        a.plot(g.cutoff, 100 * g.accuracy, "-", color=GRID, lw=1.1, zorder=2)
    a.plot(cuts, acc, "-o", color=PERF, lw=2.4, ms=5.5, zorder=4,
           markeredgecolor="white", markeredgewidth=1.1,
           label="typical accuracy across the 28 endpoints")
    a.plot(cuts, cov, "-s", color=COV, lw=2.4, ms=5, zorder=4,
           markeredgecolor="white", markeredgewidth=1.1,
           label="compounds it will answer for")
    a.set_xlabel("required similarity to a reference compound (Tanimoto)",
                 fontsize=9)
    a.set_ylabel("percent", fontsize=9)
    a.set_title("A   Accuracy for each of the 28 endpoints, in gray, as the "
                "required match tightens. The share answered keeps falling and "
                "does not level off", fontsize=9.5, color=INK, loc="left",
                pad=6)
    # The share answered keeps falling to single figures at the tightest
    # match. Starting the axis at 40% cut that curve off mid-descent and hid
    # that it never levels out, so the axis runs to zero.
    a.set_ylim(0, 103)
    a.set_yticks([0, 25, 50, 75, 100])
    a.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])
    a.legend(fontsize=8.5, frameon=False, loc="lower left")
    style(a)

    b = fig.add_axes([0.205, 0.075, 0.755, 0.475])
    b.set_facecolor("white")
    # Ranked bars, best at the top. A scatter of the two rates was tried and
    # discarded: putting them on two axes invites reading a correlation between
    # them, when they are two coordinates of one quality rather than a
    # relationship, and colouring by Matthews correlation then repeated what
    # position already said.
    op = op.sort_values("mcc").reset_index(drop=True)
    thin = op.n_toxic.fillna(0) < 20
    y = range(len(op))
    b.barh(list(y), op.mcc, height=0.62, zorder=2,
           color=[GRID if t else PERF for t in thin],
           edgecolor=[MUTED if t else "none" for t in thin], linewidth=0.9)
    b.set_yticks(list(y))
    b.set_yticklabels(op.endpoint, fontsize=8)
    for i, r in op.iterrows():
        b.text(r.mcc + 0.012, i, "%.2f" % r.mcc, va="center", ha="left",
               fontsize=7.6, color=INK, zorder=3)
        b.text(1.14, i, "catches %3.0f%%,  passes %3.0f%%"
               % (100 * r.caught, 100 * r.passed), va="center", ha="left",
               fontsize=7.4, color=MUTED, zorder=3)
    b.set_xlim(0, 1.72)
    b.set_ylim(-0.7, len(op) - 0.3)
    b.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    b.set_xlabel("Matthews correlation between the call and the measurement. "
                 "Gray bars rest on fewer than 20 measured actives.",
                 fontsize=9)
    b.set_title("B   Where the method works, at the strictest match that still "
                "answers one compound in five", fontsize=9.5, color=INK,
                loc="left", pad=6)
    style(b, grid="x")
    out = os.path.join(OUTDIR, "figure3_endpoints.png")
    fig.savefig(out, facecolor="white")
    print("wrote %s" % out)
    op.to_csv(os.path.join(OUTDIR, "endpoint_operating_points.csv"), index=False)


def main():
    res = json.load(open(os.path.join(RES, "refsim_sweeps.json")))
    fig_2003(res)
    fig_endpoints(res)


if __name__ == "__main__":
    main()
