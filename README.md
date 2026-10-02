# ogm-tandem-repeat-sizing

Estimate the size (number of repeat units) of a tandem repeat from **Bionano Optical Genome Mapping (OGM)** single-molecule alignments.

For every molecule spanning a locus of interest, the script measures the distance between the two reference labels flanking the repeat, compares it to the distance expected without expansion, and converts the difference into a number of repeat units.

It was developed for the **G4C2 hexanucleotide expansion in *C9orf72***, which is the default target, but every locus-specific parameter can be set from the command line, so the same script can be used for any other repeat expansion covered by two flanking labels.

> **Status:** research code accompanying a manuscript in preparation. Estimates are derived from label spacing and are not base-pair-resolution measurements (see [Limitations](#limitations)).

Developed at the Genomics Platform, Department of Translational Research, Institut Curie.

---

## Table of contents

- [How it works](#how-it-works)
- [Installation](#installation)
- [Usage](#usage)
- [Input files](#input-files)
- [Output](#output)
- [Using the script on another locus](#using-the-script-on-another-locus)
- [Quick test](#quick-test)
- [Limitations](#limitations)
- [Repository structure](#repository-structure)

---

## How it works

1. **Select molecules.** The `.xmap` file is filtered on the reference contig carrying the target region (`--ref-contig`).
2. **Locate the flanking labels.** In the `Alignment` field of each remaining molecule, the script looks for the consecutive reference label pair `(label-left, label-right)` and retrieves the matching query `SiteID`s. Molecules that do not contain this pair are discarded.
3. **Get positions.** The position (bp) of each query label is read from the query cmap (`_q.cmap`), which is parsed in chunks to handle large files.
4. **Compute the repeat number.** For each molecule:

```
distance_bp = Position(right label) - Position(left label)      (on the query molecule)
repeats     = (distance_bp - baseline_bp) / repeat_unit_bp
```

where `baseline_bp` is the inter-label distance in the reference without expansion and `repeat_unit_bp` is the length of the repeated motif.

If a molecule matches the label pair more than once, only the first match is kept and a warning reporting the number of such molecules is logged.

---

## Installation

Requirements: Python 3.8+ and `pandas`.

```bash
git clone https://github.com/<your-user>/ogm-tandem-repeat-sizing.git
cd ogm-tandem-repeat-sizing
pip install -r requirements.txt
```

A virtual environment (`python -m venv .venv`) or conda environment is recommended.

---

## Usage

Default target (C9orf72, G4C2):

```bash
python ogm_tandem_repeat_sizing.py \
    --xmap sample.xmap \
    --cmap sample_q.cmap \
    --output sample_repeats.csv
```

### Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--xmap` | yes | – | Alignment file (`.xmap`) |
| `--cmap` | yes | – | Query cmap (`_q.cmap`) |
| `--output` | yes | – | Output CSV path |
| `--ref-contig` | no | `9` | Reference contig ID carrying the region |
| `--label-left` | no | `6238` | Left flanking reference label |
| `--label-right` | no | `6239` | Right flanking reference label |
| `--baseline-bp` | no | `3788` | Reference inter-label distance without expansion (bp) |
| `--repeat-unit` | no | `6` | Length of the repeated motif (bp; G4C2 = 6) |

Run `python ogm_tandem_repeat_sizing.py --help` for the full list.

> **Important:** the default contig ID, label IDs and baseline are specific to the reference map used during development (**TODO: state reference genome / Bionano reference version here**). Label IDs change with the reference and the labeling enzyme, so check them against your own reference before running the script on a different setup.

---

## Input files

| File | Description |
|------|-------------|
| `.xmap` | Molecule-to-reference alignments. Must contain the columns `QryContigID`, `RefContigID` and `Alignment`. |
| `_q.cmap` | Query maps of the same alignment run. Must contain the columns `CMapId`, `SiteID` and `Position`. |

Column names are read from the `#h` header line of each file, so the column order does not matter.

---

## Output

A CSV file with one row per molecule for which the size could be estimated:

| Column | Description |
|--------|-------------|
| `QryContigID` | Identifier of the molecule (query map) |
| `distance_bp` | Distance between the two flanking labels on the molecule (bp) |
| `repeats` | Estimated number of repeat units, `(distance_bp - baseline_bp) / repeat_unit_bp` |

`repeats` is not rounded. A value near 0 corresponds to the reference-size allele; negative values indicate a distance shorter than the baseline.

---

## Using the script on another locus

To size a different repeat, provide the five target parameters:

```bash
python ogm_tandem_repeat_sizing.py \
    --xmap sample.xmap \
    --cmap sample_q.cmap \
    --output sample_locus_repeats.csv \
    --ref-contig <contig_id> \
    --label-left <label_id> \
    --label-right <label_id> \
    --baseline-bp <reference_distance_bp> \
    --repeat-unit <motif_length_bp>
```

The values can be obtained from the reference cmap used for the alignment:

- `--ref-contig`, `--label-left`, `--label-right`: the contig and the two consecutive labels (`SiteID`) that flank the repeat.
- `--baseline-bp`: the distance between these two labels in the reference (difference of their `Position` values).
- `--repeat-unit`: the motif length. If the reference allele already contains repeats, the baseline distance includes them, so the output is the **difference in repeat units relative to the reference**, not the absolute number of repeat units.

---

## Quick test

A small **synthetic** dataset is provided in [`examples/`](examples/) to check the installation:

```bash
python ogm_tandem_repeat_sizing.py \
    --xmap examples/example.xmap \
    --cmap examples/example_q.cmap \
    --output results.csv
diff results.csv examples/expected_output.csv && echo "OK"
```

Expected output:

```
QryContigID,distance_bp,repeats
101,3788.0,0.0
102,4058.0,45.0
103,4388.0,100.0
```

---

## Limitations

- Sizes are **inferred from label spacing** on optically measured molecules. Precision is bounded by OGM resolution and by stretch/sizing error, and does not reach base-pair resolution. Interpret estimates accordingly and validate against an orthogonal method where possible.
- Only molecules **spanning both flanking labels** are sized.
- When several matches of the label pair exist on a molecule, only the first is used.
- The script reports one value per molecule and does **not** assign alleles, cluster molecules or call genotypes.
- Defaults are tied to a specific reference and labeling scheme (see the note under [Arguments](#arguments)).

---

## Repository structure

```
ogm-tandem-repeat-sizing/
├── ogm_tandem_repeat_sizing.py   # main script
├── requirements.txt              # Python dependencies
├── examples/
│   ├── README.md                 # description of the test data
│   ├── example.xmap              # synthetic alignments
│   ├── example_q.cmap            # synthetic query maps
│   └── expected_output.csv       # expected result
├── CITATION.cff
├── LICENSE
├── .gitignore
└── README.md
```
