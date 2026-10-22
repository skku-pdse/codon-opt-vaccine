"""homotetramer_analysis.py — Homotetramer (mononucleotide-repeat) detection.

A *homotetramer* is defined here, following the manuscript (Table 1), as a
MAXIMAL contiguous run of four or more identical nucleotides, counted ONCE
per run — NOT once per sliding 4-nt window.

    e.g.  'AAAAAA' (6 A's) is ONE run of length 6 (counted once),
          not three overlapping 'AAAA' windows.

This is the key difference from a fixed-substring motif search (as used for
exclusion / slippery motifs): those count every overlapping occurrence, which
would over-count long runs.  Homotetramers therefore require this dedicated
run-length counter and cannot be reproduced by the substring-motif engine.

Bug fix history:
  20230926 - Fixed incorrect frame classification logic (previous version).
  Final    - Switched from overlapping 4-mer occurrence counting to maximal-run
             counting so the output matches the manuscript definition
             (each run counted once), reproducing Table 1.

Input : sequence(DNA).txt   (FASTA; RNA or DNA, any case)
Output: Homotetramer_analysis.tsv
        Columns:
          Sequence, Length,
          Number of homotetramer runs (total),
          AAAA (A-runs), TTTT (T-runs), CCCC (C-runs), GGGG (G-runs),
          Run details  (base:1-based_start:run_length, comma-separated)
"""

import re

MIN_RUN = 4  # minimum run length to be called a homotetramer


def read_fasta_file(file_path):
    """Parse a FASTA file (supports multi-line sequences).

    Sequences are upper-cased and U is converted to T so that RNA or DNA,
    in any case, are handled identically.

    Returns:
        dict: {sequence_name: sequence_string}
    """
    sequences = {}
    current = None
    with open(file_path, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                current = line[1:].strip()
                sequences[current] = ''
            elif current is not None:
                sequences[current] += ''.join(line.split()).upper().replace('U', 'T')
    return sequences


def find_homopolymer_runs(sequence, min_run=MIN_RUN):
    """Return every maximal run of >= min_run identical nucleotides.

    Args:
        sequence (str): Upper-case nucleotide sequence (A/C/G/T).
        min_run  (int): Minimum run length (default 4).

    Returns:
        list of (base, start_0based, run_length), in ascending position order.
        Runs are maximal and non-overlapping, so each is counted exactly once.
    """
    runs = []
    # ([ACGT])\1{min_run-1,}  ->  one base followed by >=(min_run-1) copies of
    # the SAME base = a maximal run of length >= min_run (greedy, non-overlapping)
    pattern = re.compile(r'([ACGT])\1{%d,}' % (min_run - 1))
    for m in pattern.finditer(sequence):
        runs.append((m.group(1), m.start(), len(m.group(0))))
    return runs


# --- Parameters ---
file_path = 'sequence(DNA).txt'
output_file = 'Homotetramer_analysis.tsv'

# --- Run analysis ---
sequences = read_fasta_file(file_path)

rows = []
for name, seq in sequences.items():
    runs = find_homopolymer_runs(seq)
    per_base = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    for base, _start, _length in runs:
        per_base[base] += 1
    details = ', '.join(f'{base}:{start + 1}:{length}' for base, start, length in runs)
    rows.append((
        name, len(seq), len(runs),
        per_base['A'], per_base['T'], per_base['C'], per_base['G'],
        details,
    ))

# --- Write output ---
with open(output_file, 'w') as fh:
    fh.write(
        'Sequence\tLength\tNumber of homotetramer runs\t'
        'AAAA (A-runs)\tTTTT (T-runs)\tCCCC (C-runs)\tGGGG (G-runs)\t'
        'Run details (base:1-based_start:length)\n'
    )
    for (name, length, total, a, t, c, g, details) in rows:
        fh.write(f'{name}\t{length}\t{total}\t{a}\t{t}\t{c}\t{g}\t{details}\n')

print(f'Results written to {output_file}.')
