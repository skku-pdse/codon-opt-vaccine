"""codon_similarity.py

Python re-implementation of the original codon_compare.pl script.

Reads a FASTA file of coding sequences and computes a pairwise codon-level
similarity matrix.  For each pair of sequences (s1, s2), similarity is defined as
the percentage of codon positions where s1 and s2 have the same codon:

    similarity(s1, s2) = 100 * (number of matching codons) / (number of codons in s1)

Both U→T conversion and upper-casing are applied automatically, so the input
FASTA may contain RNA or DNA sequences in any case.

Outputs
-------
similarity.txt
    Tab-separated similarity matrix (rows and columns correspond to sequences
    in the same order as they appear in the input FASTA; each cell is a float with 6 decimal places).
order.txt
    Sequence names in the same order as the input FASTA.

"""


def read_fasta(filename):
    """Parse a FASTA file into a dict of {name: DNA_sequence}.

    - Converts RNA to DNA (U → T).
    - Converts all characters to uppercase.
    - Strips whitespace from sequence lines.

    Args:
        filename (str): Path to the FASTA file.

    Returns:
        dict: {sequence_name (str): sequence (str)}
    """
    seq_hash = {}
    header = ''
    seq_parts = []

    with open(filename, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n').rstrip('\r')
            if line.startswith('>'):
                if header:
                    seq_hash[header] = ''.join(seq_parts)
                header = line[1:].strip()
                seq_parts = []
            else:
                processed = line.upper().replace('U', 'T').replace(' ', '')
                seq_parts.append(processed)

    if header:
        seq_hash[header] = ''.join(seq_parts)

    return seq_hash


def split_codons(sequence):
    """Split a nucleotide sequence into a list of non-overlapping codons.

    Trailing nucleotides that do not form a complete codon are discarded.

    Args:
        sequence (str): Nucleotide sequence (DNA alphabet).

    Returns:
        list: Codons as 3-character strings.
    """
    return [sequence[i:i + 3] for i in range(0, len(sequence) - len(sequence) % 3, 3)]


def compute_similarity(seq1, seq2):
    """Compute codon-level percent similarity between two sequences.

    Comparison is based on the length of seq1 (shorter comparisons
    are truncated to seq1's codon count, matching the original Perl behaviour).

    Args:
        seq1 (str): Reference nucleotide sequence.
        seq2 (str): Query nucleotide sequence.

    Returns:
        float: Percentage of codons in seq1 that are identical at the same
               position in seq2.  Returns 0.0 if seq1 has no complete codons.
    """
    codons1 = split_codons(seq1)
    codons2 = split_codons(seq2)

    if not codons1:
        return 0.0

    count = sum(1 for i, c in enumerate(codons1) if i < len(codons2) and c == codons2[i])
    return 100.0 * count / len(codons1)


def process(input_filename):
    """Run the full pairwise similarity analysis and write output files.

    Args:
        input_filename (str): Path to the input FASTA file.
    """
    seq_hash = read_fasta(input_filename)
    sequence_names = list(seq_hash.keys())

    with open('similarity.txt', 'w') as fh_sim, open('order.txt', 'w') as fh_ord:
        for name1 in sequence_names:
            fh_ord.write(f'{name1}\n')
            row_values = []
            for name2 in sequence_names:
                sim = compute_similarity(seq_hash[name1], seq_hash[name2])
                row_values.append(f'{sim:.6f}')
            fh_sim.write('\t'.join(row_values) + '\n')


# --- Parameters ---
file_path = 'sequence(DNA).txt'

# --- Run analysis ---
process(file_path)
