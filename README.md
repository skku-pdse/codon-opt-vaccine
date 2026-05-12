# JRSI Vaccine Codon Optimization — Analysis Code

Code accompanying the manuscript:

> **[Manuscript title — placeholder]**
> Je Hun Moon
> *Journal of the Royal Society Interface* 
> DOI: *[to be added upon publication]*

**Author**: Je Hun Moon
**Affiliation**: Process Design & Systems Engineering Lab, Sungkyunkwan University
**Copyright**: © 2026 Je Hun Moon
**License**: GNU General Public License v3.0 (GPL-3.0) — see [LICENSE](LICENSE)

---

## Overview

This repository contains the analysis pipeline used in the paper to evaluate
candidate codon-optimized sequences against multiple sequence-design
criteria. The pipeline is organised into five independent modules that can be
run in any order; each module reads a DNA FASTA file and writes its results to
the same folder.

| # | Module | Purpose |
|---|--------|---------|
| 01 | `01_codon_bias` | Compute reference codon usage bias (ICU & CC) from a gene set |
| 02 | `02_codon_index` | Score sequences with ICU, CC, and CAI indices |
| 03 | `03_cpg_analysis` | Detect CpG dinucleotides per reading frame |
| 04 | `04_codon_similarity` | Compute pairwise codon-level percent identity matrix |
| 05 | `05_exclusion_motif` | Locate splice / TATA-box / regulatory motifs to be avoided |

G-quadruplex propensity was assessed with the external G4Hunter tool — see
[External tool: G4Hunter](#external-tool-g4hunter) below.

---

## System requirements

- **Operating system**: Windows / macOS / Linux
- **Python**: ≥ 3.8
- **Dependencies** (see [`requirements.txt`](requirements.txt)):
  - `numpy ≥ 1.20`
  - `pandas ≥ 1.3`

```bash
pip install -r requirements.txt
```

---

## Repository layout

```
codon-opt-vaccine/
├── README.md
├── LICENSE                       # GPL-3.0
├── requirements.txt
├── 01_codon_bias/
├── 02_codon_index/
├── 03_cpg_analysis/
├── 04_codon_similarity/
└── 05_exclusion_motif/
```

Each module folder contains:
- the analysis script(s) (`*.py`),
- a copy of the required input file(s),
- and the output files produced when the script was last executed
  (these serve as reference outputs for verifying reproducibility).

---

## Running each module

All scripts are run from inside their own module folder. Input files are read
from the current working directory (filenames are hard-coded at the top of each
script — edit if you wish to use different input files).

### Module 01 — Compute codon bias
```bash
cd 01_codon_bias
python compute_icu_bias.py        # Individual codon usage (ICU)
python compute_cc_bias.py         # Codon-pair usage (CC)
```
- **Input**: `codon usage bias.txt` (TSV with columns `Gene`, `Nucleotide`),
  `trans_table.txt`, `Pair.txt`
- **Output**: `ICU_results.tsv`, `CC results.tsv`, plus per-amino-acid /
  per-codon count files

### Module 02 — Compute codon index (ICU, CC, CAI)
```bash
cd 02_codon_index
python compute_score.py
```
- **Input**: `sequence(DNA).txt` (FASTA), `trans_table.txt`,
  `count_codon.txt`, `count_aa.txt`, `count_codonp.txt`, `count_aap.txt`
- **Output**: `scores.txt` — tab-separated columns `Label`, `IC`, `CC`, `CAI`
  (14-decimal precision)
- Multi-process; uses all available CPU cores.

**Preparing the reference-table inputs from Module 01.**
The four `count_*.txt` files consumed by this module are produced by Module 01
under different names. After running `compute_icu_bias.py` and
`compute_cc_bias.py` in `01_codon_bias/`, copy and rename the output files
into `02_codon_index/` as follows:

| Source (`01_codon_bias/`) | Destination (`02_codon_index/`) |
|---------------------------|----------------------------------|
| `ICU results_Count.txt`   | `count_codon.txt`                |
| `ICU results_AA.txt`      | `count_aa.txt`                   |
| `CC results_Count.txt`    | `count_codonp.txt`               |
| `CC results_AA.txt`       | `count_aap.txt`                  |

The current repository already includes a pre-computed set of these files
generated from the reference gene list used in the manuscript; re-deriving
them is only needed when changing the reference gene set.

### Module 03 — CpG analysis
```bash
cd 03_cpg_analysis
python cpg_frame_analysis.py
```
- **Input**: `sequence(DNA).txt` (FASTA)
- **Output**: `CpG_frame_analysis.tsv` — frame-wise locations of all CpG
  dinucleotides in each input sequence.

### Module 04 — Codon-level similarity
```bash
cd 04_codon_similarity
python codon_similarity.py
```
- **Input**: `sequence(DNA).txt` (FASTA)
- **Output**: `similarity.txt` (% identity matrix), `order.txt` (row/column order)

### Module 05 — Exclusion motif detection
```bash
cd 05_exclusion_motif
python exclusion_motif_analysis.py
```
- **Input**: `sequence(DNA).txt` (FASTA)
- **Output**: `Exclusion_analysis.tsv`
- Motifs searched: `CAGG`, `AAGGTAAGT`, `AAGGTGAGT`, `CAGGTAAGT`, `CAGGTGAGT`,
  `TATAAA`, `GGCCAATCT`, `GGTCAATCT`.

---

## Output column reference

For modules 03 and 05 the TSV output columns are:

| Column | Description |
|--------|-------------|
| `Sequence` | FASTA header of the input sequence |
| `Length` | Sequence length (nt) |
| `Number of motifs (in-frame)` | Count of occurrences starting at codon position 0 (0-based index `% 3 == 0`) |
| `Motifs (in-frame)` | Comma-separated list of matched substrings |
| `Locations (in-frame)` | Comma-separated 1-based start positions (nt) |
| `Number of motifs (out-frame)` | Count of out-of-frame occurrences |
| `… (+1 frame)` / `… (+2 frame)` | Subset of out-of-frame matches at offset +1 or +2 |

For module 02 the output columns are:

| Column | Description |
|--------|-------------|
| `Label` | FASTA header of the scored sequence |
| `IC` | Individual Codon Usage score (negative; closer to 0 = better match to reference) |
| `CC` | Codon Context (codon-pair) score (negative; closer to 0 = better) |
| `CAI` | Codon Adaptation Index (0–1; higher = better) |

---

## Data

**Reference codon usage tables** (used by modules 01 and 02) were derived from
the **Human Protein Atlas** (HPA) "RNA consensus tissue gene data" — specifically
the **top 1% most highly expressed protein-coding genes in human skeletal muscle**
(ranked by descending nTPM). Mitochondrial, ribosomal, retired, and pseudogenes
were filtered out using NCBI RefSeq annotations prior to ranking. See the
accompanying manuscript and its Supplementary Material for the full reference
gene list and methodological details, including alternative reference sets
(KAZUSA, CoCoPUTs, ribosomal-only).

- **`*/sequence(DNA).txt`** — input candidate sequences in FASTA format. The
  copy in each module folder is identical and is included so each module can be
  run from its own directory without path changes.

---

## External tool: G4Hunter

G-quadruplex propensity in the candidate sequences was assessed with the
**G4Hunter** algorithm (Bedrat, Lacroix & Mergny, *Nucleic Acids Res.* 44, 1746–1759, 2016).

- Source repository: <https://github.com/AnimaTardeb/G4-hunter>
- License of G4Hunter: GNU General Public License v3.0
- Parameters used in this study: **window = 25**, **score threshold = 1.2**
  (G4Hunter defaults)

G4Hunter is **not redistributed** in this repository. To reproduce the
G-quadruplex analysis, clone G4Hunter from the original repository and run it
on the same `sequence(DNA).txt` input file used by modules 03–05.

---

## How to cite

If you use this code, please cite the accompanying manuscript:

> Moon, J. H. *[Manuscript title — placeholder]*. *J. R. Soc. Interface* (year).
> DOI: *[to be added upon publication]*.

When citing the G-quadruplex analysis, please additionally cite:

> Bedrat, A., Lacroix, L. & Mergny, J.-L. Re-evaluation of G-quadruplex
> propensity with G4Hunter. *Nucleic Acids Res.* **44**, 1746–1759 (2016).

---

## License

This project is licensed under the **GNU General Public License v3.0**.
See the [LICENSE](LICENSE) file for the full license text.

```
Copyright (C) 2026  Je Hun Moon
Process Design & Systems Engineering Lab, Sungkyunkwan University

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
```

---

## Contact

[Contact via corresponding author]
