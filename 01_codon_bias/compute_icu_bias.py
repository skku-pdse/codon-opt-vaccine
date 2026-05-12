import os
import pandas as pd
import numpy as np
from multiprocessing import Pool
import time

# Amino acid to codon mapping (standard genetic code, DNA alphabet)
# Keys: single-letter amino acid codes; '*' = stop codon
# Codon lists within each amino acid are sorted alphabetically
AminoAcids = {
    'A': ['GCA', 'GCC', 'GCG', 'GCT'],
    'C': ['TGC', 'TGT'],
    'D': ['GAC', 'GAT'],
    'E': ['GAA', 'GAG'],
    'F': ['TTC', 'TTT'],
    'G': ['GGA', 'GGC', 'GGG', 'GGT'],
    'H': ['CAC', 'CAT'],
    'I': ['ATA', 'ATC', 'ATT'],
    'K': ['AAA', 'AAG'],
    'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
    'M': ['ATG'],
    'N': ['AAC', 'AAT'],
    'P': ['CCA', 'CCC', 'CCG', 'CCT'],
    'Q': ['CAA', 'CAG'],
    'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
    'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
    '*': ['TAA', 'TAG', 'TGA'],
    'T': ['ACA', 'ACC', 'ACG', 'ACT'],
    'V': ['GTA', 'GTC', 'GTG', 'GTT'],
    'W': ['TGG'],
    'Y': ['TAC', 'TAT'],
}

# All 64 codons sorted in alphabetical order.
# NOTE: This list determines the column order in the output TSV.
#       The internal order must match the codon lists in AminoAcids above
#       (each codon must appear in exactly one amino acid group).
CodonListmanual = sorted([codon for codons in AminoAcids.values() for codon in codons])


def processGen(param):
    """Extract codons and compute ICU (Individual Codon Usage) frequencies.

    Args:
        param (str): Nucleotide sequence (DNA, length must be a multiple of 3).

    Returns:
        tuple: (CodonNormFreq, AminoCodonListFreq, CodonFreq)
            - CodonNormFreq   : dict {codon: normalised frequency within its amino acid}
            - AminoCodonListFreq : dict {amino_acid: [raw_count_per_codon, ...]}
            - CodonFreq       : dict {codon: raw_count}
    """
    codon = []               # ordered list of codons extracted from the sequence
    CodonFreq = {}           # raw count of each codon
    CodonNormFreq = {}       # normalised frequency of each codon within its amino acid
    AminoFreq = {}           # total count for each amino acid
    AminoCodonListFreq = {}  # per-codon counts grouped by amino acid
    AminoNormFreq = {}       # normalised per-codon counts grouped by amino acid

    # Step 1: split sequence into non-overlapping codons
    # Use integer division (//) because valid CDS lengths are always multiples of 3
    for i in range(len(param) // 3):
        j = i * 3
        codon.append(param[j:j + 3])

    # Step 2: count occurrences of each codon and group by amino acid
    for key in AminoAcids:
        CinA = AminoAcids[key]
        lst = []
        cumm = 0
        for k in CinA:
            cnt = codon.count(k)
            CodonFreq[k] = cnt
            lst.append(cnt)
            cumm += cnt
        AminoFreq[key] = cumm
        AminoCodonListFreq[key] = lst

    # Step 3: normalise each codon count within its amino acid
    for key in AminoCodonListFreq:
        CinA = AminoCodonListFreq[key]
        lst = []
        if int(AminoFreq[key]) != 0:
            for k in CinA:
                lst.append(round(float(k / AminoFreq[key]), 4))
        else:
            lst = CinA.copy()
        AminoNormFreq[key] = lst

    # Step 4: build a flat codon -> normalised-frequency mapping
    for key in AminoAcids:
        Clist = AminoAcids[key]
        ClistFreq = AminoNormFreq[key]
        for (c, d) in zip(Clist, ClistFreq):
            CodonNormFreq[c] = d

    return CodonNormFreq, AminoCodonListFreq, CodonFreq


def process_each_gene(nucleotide_sequence):
    """Wrapper for multiprocessing: strip newlines and call processGen."""
    param = nucleotide_sequence.replace('\n', ' ')
    return processGen(param)


def main():
    start_time = time.time()

    # Load sequence data from TSV file (columns: 'Gene', 'Nucleotide')
    # TSV format is used because Excel has a single-cell character limit of 32,767,
    # which can be exceeded by large gene sequences or full viral genomes.
    tsv_data = pd.read_csv('codon usage bias.txt', sep='\t', engine='python', encoding='utf-8')

    sequences = tsv_data['Nucleotide'].tolist()
    genes = tsv_data['Gene'].tolist()

    # Parallel processing: use all available CPU cores
    with Pool(processes=os.cpu_count()) as pool:
        results = pool.map(process_each_gene, sequences)

    # Aggregate raw codon counts and amino-acid group counts across all sequences
    aggregated_AminoCodonListFreq = {}
    aggregated_CodonFreq = {cp: 0 for cp in CodonListmanual}

    for result in results:
        _, AminoCodonListFreq, CodonFreq = result

        for amino_acid, codon_freq_list in AminoCodonListFreq.items():
            aggregated_AminoCodonListFreq[amino_acid] = (
                aggregated_AminoCodonListFreq.get(amino_acid, 0) + sum(codon_freq_list)
            )

        for codon, freq in CodonFreq.items():
            aggregated_CodonFreq[codon] += freq

    # Compute dataset-level normalised frequencies using trans_table.txt
    # trans_table.txt format: "codon:amino_acid;" repeated for all 64 codons
    normalized_CodonFreq = {}
    with open('trans_table.txt', 'r') as file:
        pair_txt = file.read().strip(';')
        pair_dict = dict(item.split(':') for item in pair_txt.split(';'))

    for codon, total_freq in aggregated_CodonFreq.items():
        amino_acid = pair_dict.get(codon)
        if amino_acid:
            total_aa_freq = aggregated_AminoCodonListFreq.get(amino_acid, 0)
            normalized_CodonFreq[codon] = total_freq / total_aa_freq if total_aa_freq else 0

    # Write intermediate count files (used for visual inspection; not required for downstream tools)
    with open('ICU results_Count.txt', 'w') as f:
        for codon, freq in aggregated_CodonFreq.items():
            f.write(f'{codon}\t{freq}\n')

    with open('ICU results_AA.txt', 'w') as f:
        for amino_acid, freq in aggregated_AminoCodonListFreq.items():
            f.write(f'{amino_acid}\t{freq}\n')

    # Build the output DataFrame: genes as rows, codons (alphabetical) as columns
    # Final two rows: 'Total' (raw sum) and 'NormTot' (normalised frequency for the whole dataset)
    excel_data = {'Gene': genes + ['Total', 'NormTot']}
    for cp in CodonListmanual:
        excel_data[cp] = (
            [r[2].get(cp, 0) for r in results]
            + [aggregated_CodonFreq[cp], round(normalized_CodonFreq.get(cp, 0), 4)]
        )

    df = pd.DataFrame(excel_data)
    df.to_csv('ICU_results.tsv', sep='\t', index=False)

    print('Execution time: ', time.time() - start_time)


if __name__ == '__main__':
    main()
