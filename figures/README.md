# The figures in the paper

Every figure in *Toxicity by Consortium Revisited: Structural Neighbors and
Cross-Family Reach* is redrawn by the scripts here, from the measured result
files in `data/`. They need matplotlib, which the base install leaves out:

```
pip install 'toxpred[figures]'      # or: pip install matplotlib
```

Run them all:

```
toxpred figures --out figures/out
```

or one at a time:

```
python figures/make_graphical_abstract.py   # Graphical abstract
python figures/make_figures_paper.py        # Figures 1 and 2
python figures/make_figure4_breadth.py      # Figure 3
python figures/make_figure6_useful.py       # Figure 4
python figures/make_figure_structures.py    # Figure S1
python figures/make_figure7_scaffold_blind.py   # Figure S2
```

| Paper | Written as | Reads |
|---|---|---|
| Graphical abstract | `figure0_graphical_abstract.png` | `refsim_sweeps.json`, `phase10_targets_breadth.json` |
| Figure 1, the consortium sweeps | `figure2_consortium_sweeps.png` | `refsim_sweeps.json` |
| Figure 2, the 28 endpoints | `figure3_endpoints.png` | `refsim_sweeps.json` |
| Figure 3, families reached | `figure4_breadth_and_targets.png` | `phase10_targets_breadth.json` |
| Figure 4, the library triage | `figure6_usefulness.png` | `phase12_usefulness.json` |
| Figure S1, twelve compounds | `figure_structures.png` | `figure_structures_provenance.csv` |
| Figure S2, the blinded screens | `figure7_scaffold_blind.png` | `phase13_scaffold_blind.json` |

The file names carry the number each figure had while the work was running,
which is the name the result files use, and the table is how they map onto the
numbering in version 2 of the paper (bioRxiv, 14 October 2026): four figures
in the main text and two in the Supporting Information.

The paper's PNGs were drawn with matplotlib 3.10.0 on Python 3.10, and that
pairing reproduces them pixel for pixel. A newer matplotlib draws the same
figures with slightly different text rasterization, and a newer RDKit lays
out the twelve structures in Figure S1 at different orientations (the paper's
were drawn with RDKit 2025.09.6). The molecules, counts and labels are the
same; those are the only differences to expect.

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
