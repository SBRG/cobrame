from __future__ import print_function, division, absolute_import

amino_acids = {
    "A": "ala__L_c",
    "R": "arg__L_c",
    "N": "asn__L_c",
    "D": "asp__L_c",
    "C": "cys__L_c",
    "E": "glu__L_c",
    "Q": "gln__L_c",
    "G": "gly_c",
    "H": "his__L_c",
    "I": "ile__L_c",
    "L": "leu__L_c",
    "K": "lys__L_c",
    "M": "met__L_c",
    "F": "phe__L_c",
    "P": "pro__L_c",
    "S": "ser__L_c",
    "T": "thr__L_c",
    "W": "trp__L_c",
    "Y": "tyr__L_c",
    "V": "val__L_c",
    "U": "cys__L_c"  # precursor to selenocysteine
}

codon_table = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S",
               "TCC": "S", "TCA": "S", "TCG": "S", "TAT": "Y", "TAC": "Y",
               "TAA": '*', "TAG": '*', "TGT": "C", "TGC": "C", "TGA": '*',
               "TGG": "W", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
               "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAT": "H",
               "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R",
               "CGA": "R", "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "I",
               "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
               "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S",
               "AGC": "S", "AGA": "R", "AGG": "R", "GTT": "V", "GTC": "V",
               "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A",
               "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
               "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

transcription_table = {"A": "utp_c", "T": "atp_c", "C": "gtp_c", "G": "ctp_c"}

base_pairs = {"A": "T", "T": "A", "G": "C", "C": "G"}


def reverse_transcribe(seq):
    return ''.join(base_pairs[i] for i in reversed(seq))


def return_frameshift_sequence(full_seq, frameshift_string):
    # Subtract 1 from start position to account for 0 indexing
    seq = ''
    for x in frameshift_string.split(','):
        left_pos, right_pos = x.split(':')
        seq += full_seq[int(left_pos) - 1: int(right_pos)]
    return seq


def extract_sequence(full_seq, left_pos, right_pos, strand):
    seq = full_seq[left_pos:right_pos]
    if strand == "+":
        return seq
    elif strand == "-":
        return reverse_transcribe(seq)
    else:
        raise ValueError("strand must be either '+' or '-'")


def get_amino_acid_sequence_from_dna(dna_seq):
    if len(dna_seq) % 3 != 0:
        raise ValueError('Gene nucleotide sequence is not a valid length')

    codons = (dna_seq[i: i + 3]
              for i in range(0, (len(dna_seq)), 3))
    amino_acid_sequence = ''.join(codon_table[i] for i in codons)
    amino_acid_sequence = amino_acid_sequence.rstrip("*")
    if "*" in amino_acid_sequence:
        amino_acid_sequence = amino_acid_sequence.replace('*', 'C')
    return amino_acid_sequence
