# The figures in the paper

Every figure in *Toxicity by Consortium Revisited: Structural Neighbors and
Cross-Family Reach* is redrawn by the scripts here, from the measured result
files in `data/`. Run them all:

```
toxpred figures --out figures/out
```

or one at a time:

```
python figures/make_figures_paper.py        # Figures 1 and 2
python figures/make_figure4_breadth.py      # Figure 3
python figures/make_figure_structures.py    # Figure 4
python figures/make_figure7_scaffold_blind.py   # Figure 5
python figures/make_figure6_useful.py       # Figure 6
```

| Paper | Written as | Reads |
|---|---|---|
| Figure 1, the consortium sweeps | `figure2_consortium_sweeps.png` | `refsim_sweeps.json` |
| Figure 2, the 28 endpoints | `figure3_endpoints.png` | `refsim_sweeps.json` |
| Figure 3, families reached | `figure4_breadth_and_targets.png` | `phase10_targets_breadth.json` |
| Figure 4, twelve compounds | `figure_structures.png` | `figure_structures_provenance.csv` |
| Figure 5, the blinded screens | `figure7_scaffold_blind.png` | `phase13_scaffold_blind.json` |
| Figure 6, the library triage | `figure6_usefulness.png` | `phase12_usefulness.json` |

The file names carry the number each figure had while the work was running,
which is the name the result files use, and the table is how they map onto the
numbering in the paper.

`data/` holds the measured values themselves, so a figure can be checked
against the number it draws without rerunning the analysis. Every one of them
came out of the run described in the paper: the reference set sweeps, the
endpoint operating points, the reverse screen against the index release of
20 September 2026, the two blinded screens, and the library triage on the held
out 30%. `endpoint_operating_points.csv` is written into the output directory
when Figure 2 is drawn, next to the figures rather than back over the input.

To draw the same figures from a fresh run of the analysis, point the scripts at
that run rather than editing them:

```
TOXPRED_FIGURE_DATA=/path/to/results TOXPRED_FIGURE_OUT=/tmp/figs \
  python figures/make_figures_paper.py
```

## Counts: check what the column counts

Never take a count from a summary column without checking what that column
counts. On 1 October 2026 a reviewer caught a measurement total that was a raw
record count: `toxric_endpoints.csv` carries `n`, the TOXRIC record count
**before** duplicate resolution, beside `unique`, the count after resolving to
one measurement per compound and endpoint. Every analysis uses `unique`. A
script summed `n`, asserted the total, and so pinned the wrong basis in place
for weeks, and the number reached a figure caption and the site.

Two habits follow:

- Where a script asserts a count, the assertion records its basis in a comment,
  and where two bases exist, record **both** so neither can be mistaken for the
  other.
- Any number taken from a README or a handoff document is recomputed from the
  data or the code before it is published. A number that has been copied once
  has usually been copied twice.
