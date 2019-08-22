import argparse
import os
import shutil

from Bio import SeqIO
from Bio.Alphabet import generic_rna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import common


def print_separator():
    """
    This simply prints a bar that is the width of the current terminal.
    """
    print('=' * shutil.get_terminal_size((80, 20)).columns, end='')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Convert codons in genetic sequence to have a lower G/C count, or to contain more frequent codons')
    parser.add_argument('--codon_table', type=str, nargs='?', default='data/codons.txt',
                        help='Space delimited file that maps codons to their respected amino acids')
    parser.add_argument('--freq', type=str, nargs='?', default='data/humanCodon_freqPerThousand.csv',
                        help='CSV file with no headers that defines the frequency of given codons in the human genome'),
    parser.add_argument('--seq', type=str, nargs='?', default='data/hngly1-mrna.fasta',
                        help='Sequence to optimize G/C content and codon frequency according to human genome'),
    parser.add_argument('--output', type=str, nargs='?', default=None,
                        help='Output of the resulting optimization in fasta format')
    parser.add_argument('--weight', type=float, nargs='?', default=0.8,
                        help='How aggressively to decrease the count of G/C nucleotides in the codons')

    args = parser.parse_args()

    seq_raw, seq_acids, codon_to_acid_with_freq_and_at_count = common.get_data(args.seq, args.codon_table, args.freq)

    better_seq = common.get_analysis(seq_acids, args.weight, codon_to_acid_with_freq_and_at_count)

    g_count, c_count, a_count, t_count = (
        better_seq.count('G'), better_seq.count('C'), better_seq.count('A'), better_seq.count('T'))

    print_separator()
    print(
        f'Optimizing for {args.weight:.2%} A/T count (decrease G/C) and {1.0 - args.weight:.2%} frequency in human genome')
    print_separator()
    print(better_seq)
    print_separator()
    print(
        f'Total Nucleotides: {len(better_seq)}. G Count: {g_count}, C Count: {c_count}, A Count: {a_count}, T Count: {t_count}')
    print_separator()
    gc_count = g_count + c_count
    print(f'G/C Count: {gc_count}')
    print(f'G/C Percentage: {gc_count / len(better_seq):.2%}')
    print_separator()
    output_path = args.output or f'{args.seq.replace(".fasta", "")}_optimized.fasta'
    # Make directories leading up if they do not exist already, writing the record will fail otherwise
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    record = SeqRecord(better_seq, id=seq_raw.id, description=seq_raw.description)
    SeqIO.write(record, output_path, 'fasta')
    print(f'Wrote output file to {output_path}')
    print_separator()
    print(f'Are both same amino acids: {seq_raw.translate().seq == better_seq.translate()}')
    print_separator()
