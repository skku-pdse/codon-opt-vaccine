def classify_motif_frame(sequence, motif, position):
    """Classify the reading frame of a motif occurrence.

    Args:
        sequence (str): Full nucleotide sequence.
        motif    (str): Motif string to classify.
        position (int): 0-based index of the motif start in the sequence.

    Returns:
        str: 'in-frame' if position % 3 == 0, else 'out-frame'.

    Notes:
        Python indices are 0-based, so codon boundaries are at positions
        0, 3, 6, ... (i.e., 3k for k = 0, 1, 2, ...).
        An in-frame motif starts exactly at a codon boundary (index % 3 == 0).
        Reported 1-based locations are obtained by adding 1 to the 0-based index.
    """
    index = sequence.find(motif, position)
    if index % 3 == 0:
        return 'in-frame'
    return 'out-frame'


def count_exclusion_motifs(sequence, exclusion_motifs):
    """Search for all motif occurrences and classify them by reading frame.

    Overlapping occurrences are counted individually (search advances by 1 each step).

    Args:
        sequence        (str):  Nucleotide sequence.
        exclusion_motifs (list): List of motif strings to search for.

    Returns:
        tuple of 12 elements:
            (count_in_frame, motifs_in_frame, locations_in_frame,
             count_out_frame, motifs_out_frame, locations_out_frame,
             count_plus1_frame, motifs_plus1_frame, locations_plus1_frame,
             count_plus2_frame, motifs_plus2_frame, locations_plus2_frame)
        Locations are reported as 1-based positions.
    """
    count_in_frame_motifs = 0
    count_out_frame_motifs = 0
    count_plus1_frame_motifs = 0
    count_plus2_frame_motifs = 0
    motifs_in_frame = []
    motifs_out_frame = []
    motifs_plus1_frame = []
    motifs_plus2_frame = []
    locations_in_frame = []
    locations_out_frame = []
    locations_plus1_frame = []
    locations_plus2_frame = []

    for motif in exclusion_motifs:
        pos = 0
        while pos < len(sequence):
            next_pos = sequence.find(motif, pos)
            if next_pos == -1:
                break
            motif_frame = classify_motif_frame(sequence, motif, next_pos)
            if 'in-frame' in motif_frame:
                count_in_frame_motifs += 1
                motifs_in_frame.append(sequence[next_pos:next_pos + len(motif)])
                locations_in_frame.append(next_pos + 1)  # convert to 1-based
            elif 'out-frame' in motif_frame:
                count_out_frame_motifs += 1
                motifs_out_frame.append(sequence[next_pos:next_pos + len(motif)])
                locations_out_frame.append(next_pos + 1)
                if next_pos % 3 == 1:
                    count_plus1_frame_motifs += 1
                    motifs_plus1_frame.append(sequence[next_pos:next_pos + len(motif)])
                    locations_plus1_frame.append(next_pos + 1)
                elif next_pos % 3 == 2:
                    count_plus2_frame_motifs += 1
                    motifs_plus2_frame.append(sequence[next_pos:next_pos + len(motif)])
                    locations_plus2_frame.append(next_pos + 1)
            pos = next_pos + 1

    return (
        count_in_frame_motifs, motifs_in_frame, locations_in_frame,
        count_out_frame_motifs, motifs_out_frame, locations_out_frame,
        count_plus1_frame_motifs, motifs_plus1_frame, locations_plus1_frame,
        count_plus2_frame_motifs, motifs_plus2_frame, locations_plus2_frame
    )


def get_sequence_length(sequence):
    return len(sequence)


def read_fasta_file(file_path):
    """Parse a FASTA file (supports multi-line sequences).

    Args:
        file_path (str): Path to the FASTA file.

    Returns:
        dict: {sequence_name: sequence_string}
    """
    sequences = {}
    current_sequence = ''
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                current_sequence = line[1:].strip()
                sequences[current_sequence] = ''
            else:
                sequences[current_sequence] += line.strip()
    return sequences


# --- Parameters ---
file_path = 'sequence(DNA).txt'

# CpG dinucleotide: the primary motif of interest in this analysis.
# 'CG' spans a codon boundary when it appears at the last nucleotide of one codon
# and the first of the next (inter-codon CpG), or fully within a codon (intra-codon CpG).
exclusion_motifs = ['CG']

# --- Run analysis ---
sequences = read_fasta_file(file_path)

results = []
for sequence_name, sequence in sequences.items():
    (
        count_in_frame_motifs, motifs_in_frame, locations_in_frame,
        count_out_frame_motifs, motifs_out_frame, locations_out_frame,
        count_plus1_frame_motifs, motifs_plus1_frame, locations_plus1_frame,
        count_plus2_frame_motifs, motifs_plus2_frame, locations_plus2_frame
    ) = count_exclusion_motifs(sequence, exclusion_motifs)
    sequence_length = get_sequence_length(sequence)
    results.append((
        sequence_name, sequence_length,
        count_in_frame_motifs, motifs_in_frame, locations_in_frame,
        count_out_frame_motifs, motifs_out_frame, locations_out_frame,
        count_plus1_frame_motifs, motifs_plus1_frame, locations_plus1_frame,
        count_plus2_frame_motifs, motifs_plus2_frame, locations_plus2_frame
    ))

# --- Write output ---
output_file = 'CpG_frame_analysis.tsv'
with open(output_file, 'w') as file:
    file.write(
        'Sequence\tLength\t'
        'Number of motifs (in-frame)\tMotifs (in-frame)\tLocations (in-frame)\t'
        'Number of motifs (out-frame)\tMotifs (out-frame)\tLocations (out-frame)\t'
        'Number of motifs (+1 frame)\tMotifs (+1 frame)\tLocations (+1 frame)\t'
        'Number of motifs (+2 frame)\tMotifs (+2 frame)\tLocations (+2 frame)\n'
    )
    for result in results:
        (
            sequence_name, sequence_length,
            count_in_frame_motifs, motifs_in_frame, locations_in_frame,
            count_out_frame_motifs, motifs_out_frame, locations_out_frame,
            count_plus1_frame_motifs, motifs_plus1_frame, locations_plus1_frame,
            count_plus2_frame_motifs, motifs_plus2_frame, locations_plus2_frame
        ) = result
        file.write(
            f'{sequence_name}\t{sequence_length}\t'
            f'{count_in_frame_motifs}\t{", ".join(motifs_in_frame)}\t{", ".join(str(l) for l in locations_in_frame)}\t'
            f'{count_out_frame_motifs}\t{", ".join(motifs_out_frame)}\t{", ".join(str(l) for l in locations_out_frame)}\t'
            f'{count_plus1_frame_motifs}\t{", ".join(motifs_plus1_frame)}\t{", ".join(str(l) for l in locations_plus1_frame)}\t'
            f'{count_plus2_frame_motifs}\t{", ".join(motifs_plus2_frame)}\t{", ".join(str(l) for l in locations_plus2_frame)}\n'
        )

print(f'Results written to {output_file}.')
