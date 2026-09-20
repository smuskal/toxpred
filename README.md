# toxpred

**A compound's toxicity can be gauged from the toxicities of the compounds it resembles.**

That idea is from a 2003 paper. This package rebuilds it on data anyone can
download, extends it to 28 further toxicity endpoints, and adds a second signal
the original could not compute: how many distinct protein families a compound's
chemistry reaches in the crystallographic record.

Nothing here needs a protein structure, a docking run or a conformer ensemble.

---

## What it does

```
$ toxpred score --smiles "CC(=O)Oc1ccccc1C(=O)O"

ACUTE ORAL TOXICITY, red at or under 50 mg/kg, green over 2000
  CC(=O)Oc1ccccc1C(=O)O: GREEN, 12 voters agreeing 58%, nearest 1.000

FURTHER ENDPOINTS
  endpoint                   call     agreement   voters
  CYP3A4                     clear         100%       11
  CYP2D6                     clear         100%        4
  Ames Mutagenicity          clear         100%        4
  Hepatotoxicity             TOXIC          67%        3
  Clinical toxicity          no call                   2
```

**`no call` is a real answer.** The method predicts only where the reference set
holds a compound similar enough to the query. For anything unlike what has been
measured, it declines instead of guessing, and a tie is not a majority.

---

## Install

Self-contained. Creates its own environment and touches nothing else.

```bash
git clone https://github.com/smuskal/toxpred.git
cd toxpred
conda env create -f environment.yml
conda activate toxpred
```

Then pull the public reference data once. It lands in `~/.toxpred`, or wherever
you point `TOXPRED_HOME`, never into your project.

```bash
toxpred fetch
```

That downloads roughly 90 MB: the curated rat oral acute toxicity set
distributed with OPERA, and the thirty standardized endpoint tables from TOXRIC.
Both are public.

---

## The trade you are always making

Demanding a closer match raises accuracy and lowers the share of compounds you
get an answer for. There is no setting that escapes it, so both numbers are
always reported.

![Accuracy against coverage](docs/consortium_sweeps.png)

Growing the reference set is the one move that does not cost anything. Panel C:
from 171 to 8,599 reference compounds, accuracy holds between 85 and 92 percent
while the share of compounds answered rises from 3 to 47 percent. **A bigger
shared reference set buys reach, not accuracy**, which is why pooling is worth
something to everyone who contributes.

Turn the dials yourself:

```bash
toxpred score --input molecules.smi --cutoff 0.65 --min-members 3
```

---

## The same machinery on 28 more endpoints

Five cytochrome P450 isoforms, cardiotoxicity at four thresholds, Ames
mutagenicity, carcinogenicity, hepatotoxicity, twelve nuclear receptor and
stress response assays, eye irritation and corrosion, respiratory,
developmental, reproductive and clinical toxicity.

```bash
toxpred endpoints                       # list them
toxpred score --smiles "..." --endpoints "CYP450_CYP3A4" "Hepatotoxicity_Hepatotoxicity"
```

![Twenty-eight endpoints](docs/endpoints.png)

Two of them, the estrogen receptor assay and cardiotoxicity at 30 micromolar,
are so imbalanced that inheriting labels cannot beat always answering the same
way. They are shown that way rather than left out.

---

## Cross-family reach, the second signal

How many distinct protein families does a compound's chemistry reach among the
ligands solved in protein structures? Compounds that reach widely are more often
toxic.

![Reaching more families](docs/breadth.png)

A compound whose matches reach no family inhibits CYP3A4 3 percent of the time.
One reaching sixteen or more inhibits it 51 percent of the time. The same
ordering holds for the other cytochromes, for nuclear receptor binding and for
clinical toxicity.

This is measured with no access to how many assay panels a compound has been
through, which is the confound that makes measured promiscuity hard to read.

### This part runs locally on PharmCast

The published index stores a pharmacophore fingerprint per ligand, so a query
has to be fingerprinted the same way to be comparable. PharmCast does that, it
is open, and it runs on your own machine. CPU only; the network is small enough
that a GPU buys nothing.

- **https://github.com/smuskal/pharmcast**
- **https://pharmcast.ai**

```bash
pip install git+https://github.com/smuskal/pharmcast.git
toxpred score --smiles "CC(=O)Oc1ccccc1C(=O)O" --reach
```

That is the whole setup. On first use toxpred pulls two more things and verifies
both, then works offline:

| what | from | verified against |
|---|---|---|
| the PharmCast checkpoint | `pharmcast.ai/models` | its published `SHA256SUMS` |
| the Reverse Screen index | `reversescreen.ai` download API | the per-file SHA-256 in its manifest |

```
CROSS-FAMILY REACH, index version 2026-09-20
  query                                          matches   targets  best sim
  CC(=O)Oc1ccccc1C(=O)O                               10        36     0.707
```

Every fingerprint operation is PharmCast's own: `read_pfp` for the index,
`PharmCast.words_batch` for the queries, `pharmtan_matrix` for the comparison.
None of the packing or bit ordering is reimplemented here, so the query and the
index cannot drift apart.

Group the targets into families your own way with `--family-map`, a two column
accession and family file. Without one, reach is reported as distinct targets,
which needs nothing but the index.

Everything else in this package works without PharmCast installed, and says so
rather than failing.

#### Verified from a clean machine

This path was tested from an empty virtual environment, installing only from
the public repositories, and it reproduces the numbers above exactly:

```bash
python -m venv env && source env/bin/activate
pip install git+https://github.com/smuskal/pharmcast.git
pip install git+https://github.com/smuskal/toxpred.git
toxpred score --smiles "CC(=O)Oc1ccccc1C(=O)O" --reach
```

---

## What it is good for, and what it is not

Used to triage a virtual library, the score sets aside a tenth of it. Out of
every 100 compounds in that tenth, far more carry the liability than in the
library it came from: **6 times more androgen receptor ligands, 3 times more
CYP2D6 inhibitors, 3 times more clinically toxic compounds**.

It earns its place on rare outcomes. Where half a library already carries an
endpoint there is little left to sort, and hepatotoxicity moves only from 50 in
100 to 75 in 100. A filter cannot enrich what is already everywhere.

---

## Data sources, all public

| source | what | where |
|---|---|---|
| CATMoS | rat oral acute toxicity, curated by NICEATM and US EPA | distributed with [OPERA](https://github.com/NIEHS/OPERA) |
| TOXRIC | 30 standardized toxicity endpoints | [toxric.bioinforai.tech](https://toxric.bioinforai.tech/) |
| Reverse Screen index | ligands solved in protein structures, with their targets | [reversescreen.ai](https://reversescreen.ai) |
| PharmCast | the fingerprint the index is built on | [pharmcast.ai](https://pharmcast.ai) |

Only experimental values are used as reference values. The CATMoS file also
carries consensus model predictions for compounds that were never measured, and
those are filtered out rather than inherited from.

---

## Citing

The method:

> Muskal SM, Jha SK, Kishore MP, Tyagi P. A simple and readily integratable
> approach to toxicity prediction. *J Chem Inf Comput Sci* 2003;43:1673-1678.

The data and the retrieval index carry their own citations, listed above.

---

## License

MIT. See [LICENSE](LICENSE).
