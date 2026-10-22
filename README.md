# JRSI Vaccine Codon Optimization — Analysis Code

Code accompanying the manuscript:

> **A multi-objective sequence engineering framework for rational genetic vaccine design**
> Je Hun Moon, Seo-Young Park, Kok Siong Ang, Yoon Ki Kim, Dong-Yup Lee
> *Journal of the Royal Society Interface*
> DOI: *[to be added upon publication]*

**Corresponding author**: Dong-Yup Lee (dongyuplee@skku.edu)
**Affiliation**: School of Chemical Engineering, Sungkyunkwan University, Suwon, Republic of Korea
**License**: GNU General Public License v3.0 (GPL-3.0) — see [LICENSE](LICENSE)

---

## Overview

This repository contains the descriptive analysis pipeline used to evaluate the
wild-type SARS-CoV-2 spike coding sequence, three deployed COVID-19 vaccine
sequences, and four academic candidate sequences against multiple
sequence-design criteria. The pipeline is organised into seven independent
modules; each module reads a DNA FASTA file and writes its results to the same
folder.

| # | Module | Purpose | Manuscript output |
|---|--------|---------|-------------------|
| 01 | `01_codon_bias` | Build the reference codon-usage bias (ICU & CC) from a gene set | Methods §4.2 |
| 02 | `02_codon_index` | Score sequences with ICU, CC, and CAI (COOL-compatible) | figure 3A, table S4 |
| 03 | `03_cpg_analysis` | Count CpG dinucleotides per reading frame | table 1 (CpG); CpG positions plotted in figure 1C |
| 04 | `04_codon_similarity` | Pairwise codon-level percent-identity matrix | figure 1A |
| 05 | `05_exclusion_motif` | Splice / TATA / CAAT regulatory motifs to be avoided | table 1 (Exclusion Motif) |
| 06 | `06_homotetramer` | Homotetrameric runs (maximal runs ≥ 4 nt) | table 1 (Homotetramers) |
| 07 | `07_slippery_motif` | Ribosomal slippery sequences (−1 PRF motifs) | table 1 (Slippery Motif) |

G-quadruplex propensity was assessed with the external **G4Hunter web
application** — see
[External tool: G4Hunter](#external-tool-g4hunter).

> **Sequence labels.** The eight sequences are labelled exactly as in the
> manuscript: `WT`, `BNT162b2`, `mRNA-1273`, `AZD1222`, `LD1`, `LD2`, `OG`, `mRid`.

---

## System requirements

- **Operating system**: Windows / macOS / Linux
- **Python**: ≥ 3.8 (tested on 3.9)
- **Dependencies** (see [`requirements.txt`](requirements.txt)): `numpy ≥ 1.20`, `pandas ≥ 1.3`

```bash
pip install -r requirements.txt
```

Modules 03–07 use only the Python standard library.

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
├── 05_exclusion_motif/
├── 06_homotetramer/
└── 07_slippery_motif/
```

Each module folder contains the analysis script(s) (`*.py`), a copy of the
required input file(s), and the reference output files produced when the script
was last executed (for reproducibility verification).

---

## Running each module

All scripts are run from inside their own module folder; input filenames are
hard-coded at the top of each script.

### Module 01 — Build codon bias
```bash
cd 01_codon_bias
python compute_icu_bias.py        # Individual codon usage (ICU)
python compute_cc_bias.py         # Codon-pair usage (CC)
```
- **Input**: `codon usage bias.txt` (TSV: `Gene`, `Nucleotide`), `trans_table.txt`, `Pair.txt`
- **Output**: `ICU_results.tsv`, `CC results.tsv`, plus per-amino-acid / per-codon count files
- The CC reference excludes codon-pairs whose **first** codon is a stop codon,
  matching the COOL convention used by module 02 (a stop codon has no meaningful
  "next codon" within an ORF). `Codonpairlist` therefore holds 3,904 pairs
  (64 × 64 = 4,096 minus the 3 × 64 = 192 stop-initiated pairs).
- The TSV input format is used because a single Excel cell is limited to 32,767
  characters, which long CDS exceed.

**Preparing the reference-table inputs for Module 02.**
The four `count_*.txt` files consumed by module 02 are produced by module 01
under different names. After running both scripts in `01_codon_bias/`, copy and
rename the outputs into `02_codon_index/`:

| Source (`01_codon_bias/`) | Destination (`02_codon_index/`) |
|---------------------------|----------------------------------|
| `ICU results_Count.txt`   | `count_codon.txt`                |
| `ICU results_AA.txt`      | `count_aa.txt`                   |
| `CC results_Count.txt`    | `count_codonp.txt`               |
| `CC results_AA.txt`       | `count_aap.txt`                  |

The repository already includes a pre-computed set generated from the reference
gene list used in the manuscript; re-deriving them is only needed when changing
the reference gene set.

### Module 02 — Codon index (ICU, CC, CAI)
```bash
cd 02_codon_index
python compute_score.py
```
- **Input**: `sequence(DNA).txt` (FASTA), `trans_table.txt`, `count_codon.txt`, `count_codonp.txt`
- **Also read**: `count_aa.txt`, `count_aap.txt` — consistency check only. A mismatch
  with the totals re-derived from `count_codon.txt` / `count_codonp.txt` is printed
  as a `[WARN]` line.
- **Output**: `scores.txt` — tab-separated `Label`, `IC`, `CC`, `CAI` (14-decimal precision)
- Implements the COOL codon indices (Chung & Lee, *BMC Syst. Biol.* 2012, 6:134);
  the bundled count files are the module-01 skeletal-muscle reference. Multi-process.
- CAI is computed as the geometric mean over all codons of the sequence.
- Sequences ending in two or more consecutive stop codons are normalised to a
  single terminal stop before scoring.

### Module 03 — CpG analysis
```bash
cd 03_cpg_analysis
python cpg_frame_analysis.py
```
- **Input**: `sequence(DNA).txt`
- **Output**: `CpG_frame_analysis.tsv`
- Intra-codon CpGs are those wholly inside one codon (CGN + NCG; the in-frame and
  +1-frame classes); inter-codon CpGs span a codon boundary (the +2-frame class).

### Module 04 — Codon-level similarity
```bash
cd 04_codon_similarity
python codon_similarity.py
```
- **Input**: `sequence(DNA).txt`
- **Output**: `similarity.txt` (% identity matrix), `order.txt` (row/column order)
- When lengths differ, similarity is normalised by the length of the **left-hand
  (row)** sequence, so the matrix is intentionally asymmetric (Methods §4.1).

### Module 05 — Exclusion motif detection
```bash
cd 05_exclusion_motif
python exclusion_motif_analysis.py
```
- **Input**: `sequence(DNA).txt`
- **Output**: `Exclusion_analysis.tsv`
- Motifs: `CAGG`, `AAGGTAAGT`, `AAGGTGAGT`, `CAGGTAAGT`, `CAGGTGAGT`, `TATAAA`, `GGCCAATCT`, `GGTCAATCT`
  (i.e. the 3′ acceptor `CAGG`, the 5′ donor `MAGGTRAGT`, the TATA box, and the CAAT box `GGYCAATCT`).
- The **table 1 count for a motif = (in-frame) + (out-frame)** occurrences of that motif.

### Module 06 — Homotetramer detection
```bash
cd 06_homotetramer
python homotetramer_analysis.py
```
- **Input**: `sequence(DNA).txt`
- **Output**: `Homotetramer_analysis.tsv`
- A homotetramer is a maximal contiguous run of ≥ 4 identical nucleotides, counted
  once per run — the definition used in table 1.

### Module 07 — Slippery motif detection
```bash
cd 07_slippery_motif
python slippery_motif_analysis.py
```
- **Input**: `sequence(DNA).txt`
- **Output**: `Slippery_motif_analysis.tsv`
- The eleven motifs are taken from table 1 of Caliskan, Peske & Rodnina,
  *Trends Biochem. Sci.* **40**, 265–274 (2015), transcribed to the DNA alphabet.
  The same substring-search engine as module 05 is used with this motif list.
- **Frame convention**: a 0-based start index with `% 3 == 2` places the motif in
  the −1 frame (the functional configuration for −1 PRF, because the leading X of
  the canonical X_XXY_YYZ heptamer is the last base of the preceding codon). This
  is the `+2 frame` column of the output. Table 1 reports these hits: `TTTAAAC`
  in WT (nt 2874) and `GGGAAAC` in LD2 (nt 549).

---

## Output column reference

**Modules 03, 05, 07** (substring / dinucleotide search) — one row per sequence:

| Column | Description |
|--------|-------------|
| `Sequence`, `Length` | FASTA header; sequence length (nt) |
| `Number of motifs (in-frame)` + list + 1-based locations | occurrences starting at a codon boundary (0-based index `% 3 == 0`) |
| `Number of motifs (out-frame)` + list + locations | occurrences not on a codon boundary |
| `… (+1 frame)` / `… (+2 frame)` | out-of-frame occurrences at offset +1 / +2 |

**Module 06** (maximal-run homotetramers) — one row per sequence:

| Column | Description |
|--------|-------------|
| `Sequence`, `Length` | FASTA header; sequence length (nt) |
| `Number of homotetramer runs` | total maximal runs ≥ 4 nt |
| `AAAA / TTTT / CCCC / GGGG (…-runs)` | run count per nucleotide (**these are the table 1 values**) |
| `Run details` | `base:1-based_start:run_length`, comma-separated |

**Module 02**:

| Column | Description |
|--------|-------------|
| `Label` | FASTA header of the scored sequence |
| `IC` | Individual Codon Usage score (negative; closer to 0 = better) |
| `CC` | Codon Context (codon-pair) score (negative; closer to 0 = better) |
| `CAI` | Codon Adaptation Index (0–1; higher = better) |

All modules normalise their FASTA input identically: upper-case, U → T, and all
whitespace removed.

---

## Correspondence to the manuscript

Running modules 02–07 on the bundled `sequence(DNA).txt` reproduces, exactly:

- **Table 1** — length, CpG counts, exclusion-motif counts, slippery motifs and
  homotetramer (maximal-run) counts for all eight sequences;
- **Table S4, skeletal-muscle block** — the ICU/CC/CAI values in `scores.txt`;
- **Table S5** — the CAI/ICU/CC columns, rounded, for the three clinical vaccines;
- **Figure 1A** — the similarity matrix in `similarity.txt`;
- **Figure 3A** — the ICU/CC coordinates of the eight benchmark sequences.

GC content (table 1) is `(#G + #C) / length`; it is a direct property of each
sequence and is not produced by a separate script.

**Table S2** (per-sequence codon composition) is produced by module 01: run
`compute_icu_bias.py` with the eight sequences supplied in the module-01 TSV
input format (`Gene`, `Nucleotide`) instead of the reference gene list. The
`ICU_results.tsv` gene × codon matrix is table S2 transposed.

---

## Data

**Reference codon-usage tables** (modules 01–02) were derived from the
**Human Protein Atlas** "RNA consensus tissue gene data", specifically the
**top 1% most highly expressed protein-coding genes in human skeletal muscle**
(descending nTPM), after removing mitochondrial, ribosomal, retired, and
pseudogenes using NCBI RefSeq annotations — 160 genes, 108,923 codons. See the
manuscript Methods §4.2 for the alternative reference sets (KAZUSA, CoCoPUTs,
ribosomal-protein, genome-wide HPA) reported in table S4; all of them are built
with the same module-01 procedure from the corresponding gene list.

`*/sequence(DNA).txt` — the eight candidate/benchmark sequences in FASTA format.
The copy in each module folder is identical.

> **Scope.** This repository reproduces the *descriptive* analyses (table 1,
> tables S2/S4/S5; figures 1A, 3A). The multi-objective / CAI optimisation of
> figure 3 was performed with a separately developed codon-optimisation method
> (COOL) and is not part of this package; all optimisation objectives, parameter
> settings, selected candidate sequences and objective-score tables are provided
> in the manuscript and its electronic supplementary material.

---

## External tool: G4Hunter

G-quadruplex propensity was assessed with the **G4Hunter web application**
(Brázda et al., *Bioinformatics* 35, 3493–3495, 2019), which implements the
G4Hunter algorithm of Bedrat, Lacroix & Mergny (*Nucleic Acids Res.* 44,
1746–1759, 2016).

- Server: <https://bioinformatics.ibp.cz/#/analyse/g4hunter>
- Parameters: **window = 25**, **score threshold = 1.2**
- Input: the eight sequences in `02_codon_index/sequence(DNA).txt`

The counts given in the manuscript are the server's `Quadruplexes found` values.
Note that the standalone Python implementation
(<https://github.com/AnimaTardeb/G4-hunter>) applies a different window-merging
rule and returns slightly different counts on the same input.

---

## How to cite

> Moon, J. H., Park, S.-Y., Ang, K. S., Kim, Y. K. & Lee, D.-Y.
> *A multi-objective sequence engineering framework for rational genetic vaccine
> design.* *J. R. Soc. Interface* (year). DOI: *[to be added]*.

When citing the G-quadruplex analysis, additionally cite Brázda et al. (2019)
and Bedrat et al. (2016).

---

## License

GNU General Public License v3.0. See [LICENSE](LICENSE).

```
Copyright (C) 2026  Je Hun Moon, Dong-Yup Lee
Process Design & Systems Engineering Lab, Sungkyunkwan University
```

---

## Changelog

***

- **v1.0.2 (2026-08-17)**
  - **Added `06_homotetramer/` and `07_slippery_motif/`** (scripts, inputs and
    reference outputs). These generate the `Homotetramers` and `Slippery Motif`
    columns of table 1, which were not reproducible from v1.0.1.
  - Documentation updated: author list and citation block completed,
    manuscript-output mapping added, supplementary-table numbering aligned with
    the published order, and the G4Hunter section corrected to the web
    application that was actually used.
  - Input handling unified across modules (upper-case, U → T, whitespace
    removed); reference-table consistency checks added.

  *Reference outputs for modules 01–05 are byte-identical to v1.0.1.*

***

- **v1.0.1 (2026-05-12)** — initial version-of-record deposit (modules 01–05).
