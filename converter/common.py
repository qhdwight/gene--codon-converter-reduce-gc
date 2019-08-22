import pandas as pd
from Bio import SeqIO
from Bio.Alphabet import generic_rna
from Bio.Seq import Seq


def read_codon_table(codon_table_source):
    return pd.read_csv(codon_table_source, names=['codon', 'amino_acid'], sep=' ', header=None)


def read_codon_freq(codon_freq_source):
    return pd.read_csv(codon_freq_source, usecols=[0, 2])


def read_sequence(seq_file_source):
    return SeqIO.read(seq_file_source, 'fasta')


def get_data(seq_file, codon_table_file='data/codons.txt', freq_file='data/humanCodon_freqPerThousand.csv'):
    # We assume that the csv has no header and is separated by spaces
    codons_to_acid = read_codon_table(codon_table_file)
    # Use only first and third columns which are codons and frequency count in human genome
    codon_freq = read_codon_freq(freq_file)
    seq_raw = read_sequence(seq_file)
    # Get the amino acids from the raw nucleotides so we can reverse search for better codons
    seq_acids = seq_raw.translate()
    # Combine the codons to acid table with the frequency count in the human genome for each codon
    # This allows us to cross-reference counts, codons, and the amino acids they code for
    codon_to_acid_with_freq_and_at_count = pd.merge(codons_to_acid, codon_freq)
    # Add the A and T count to each codon, our target is to maximize them instead of having G and C's
    codon_to_acid_with_freq_and_at_count['at_count'] = codon_to_acid_with_freq_and_at_count.codon.apply(
        lambda codon: codon.count('A') + codon.count('T'))
    return seq_raw, seq_acids, codon_to_acid_with_freq_and_at_count


def get_scores(df, label):
    """
    Annotate a given column in a dataframe with a score.
    1.0 will be the highest value, while 0.0 will be the lowest value. It is essentially a re-map
    :param df: The dataframe to operate on
    :param label: The name of the column
    :return: None, this function modifies the given dataframe.
    """
    max_number = df[label].max()
    min_number = df[label].min()
    range_number = max_number - min_number
    df[f'{label}_score'] = 1.0 if range_number == 0 else df[label].apply(
        lambda val: (val - min_number) / range_number)


def get_analysis(seq_acids, weight, codon_to_acid_with_freq_and_at_count):
    # Clamp value between 0.0 to 1.0
    weight = max(0.0, min(weight, 1.0))
    better_codons = ""
    # Iterate over each acid, find the best codon sequence, and add it to our new gene builder
    for acid in seq_acids:
        if acid == '*':
            # biopython represents stop as asterisk while in our codons table it is represented as "Stop"
            acid = "Stop"
        # Find codons that can represent this amino acid. We must drop NaN because by default pandas leaves them in.
        potential_codons = codon_to_acid_with_freq_and_at_count[
            codon_to_acid_with_freq_and_at_count['amino_acid'] == acid].copy()
        get_scores(potential_codons, 'number')
        get_scores(potential_codons, 'at_count')
        # Compute the weighted average between count frequency and lowest G/C content
        potential_codons['overall_score'] = potential_codons.apply(
            lambda row: (row['number_score'] * (1.0 - weight) + row['at_count_score'] * weight), axis=1)
        idx = potential_codons['overall_score'].idxmax()
        # Search back and find what codons corresponds to this index with maximum overall score
        better_codons += codon_to_acid_with_freq_and_at_count.loc[idx].codon
    better_seq = Seq(better_codons, generic_rna)
    return better_seq
