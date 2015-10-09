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
    "V": "val__L_c"
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


def extract_sequence(full_seq, left_pos, right_pos, strand):
    seq = full_seq[left_pos:right_pos]
    if strand == "+":
        return seq
    elif strand == "-":
        return reverse_transcribe(seq)
    else:
        raise ValueError("strand must be either '+' or '-'")


tRNA_to_codon = {'b2691': ['CGU', 'CGC', 'CGA'],
                 'b2692': ['CGU', 'CGC', 'CGA'],
                 'b2693': ['CGU', 'CGC', 'CGA'],
                 'b2694': ['CGU', 'CGC', 'CGA'], 'b2695': ['AGC', 'AGU'],
                 'b2590': ['GAA', 'GAG'], 'b0244': ['ACG'],
                 'b0883': ['UCC', 'UCU'], 'b1231': ['UAC', 'UAU'],
                 'b1230': ['UAC', 'UAU'], 'b4270': ['UUG'],
                 'b3853': ['GCA', 'GCG', 'GCU'], 'b3852': ['AUC', 'AUU'],
                 'b4164': ['GGC', 'GGU'], 'b4165': ['GGC', 'GGU'],
                 'b4370': ['CUG'], 'b4163': ['GGC', 'GGU'],
                 'b2189': ['CCC', 'CCU'], 'b1032': ['UCC', 'UCU'],
                 'b2348': ['AGG'], 'b3277': ['AUC', 'AUU'],
                 'b3276': ['GCA', 'GCG', 'GCU'], 'b3273': ['ACC', 'ACU'],
                 'b1986': ['AAC', 'AAU'], 'b0666': ['AUG'], 'b0665': ['CAA'],
                 'b0664': ['CAA'], 'b3969': ['GAA', 'GAG'], 'b4368': ['CUG'],
                 'b4369': ['CUG'], 'b2864': ['GGG'],
                 'b0971': ['UCU', 'UCA', 'UCG'], 'b0206': ['GAC', 'GAU'],
                 'b3757': ['GAA', 'GAG'], 'b2404': ['AAA', 'AAG'],
                 'b2401': ['GUA', 'GUG', 'GUU'],
                 'b2403': ['GUA', 'GUG', 'GUU'], 'b2396': ['GCC'],
                 'b3978': ['GGA'], 'b3979': ['ACC', 'ACU'], 'b2397': ['GCC'],
                 'b2402': ['GUA', 'GUG', 'GUU'],
                 'b3976': ['ACU', 'ACA', 'ACG'], 'b3977': ['UAC', 'UAU'],
                 'b3171': ['START'], 'b3174': ['CUC', 'CUU'], 'b1975': ['UCG'],
                 'b1977': ['AAC', 'AAU'], 'b3761': ['UGG'],
                 'b3760': ['GAC', 'GAU'], 'b4134': ['UUC', 'UUU'],
                 'b3545': ['CCG'], 'b0668': ['CAG'], 'b1984': ['AAC', 'AAU'],
                 'b0203': ['GCA', 'GCG', 'GCU'], 'b0202': ['AUC', 'AUU'],
                 'b1989': ['AAC', 'AAU'], 'b2967': ['UUC', 'UUU'],
                 'b3069': ['AUA'], 'b1909': ['UUG', 'UUA'], 'b0670': ['CAG'],
                 'b0672': ['CUG', 'CUA'], 'b0673': ['AUG'], 'b2816': ['START'],
                 'b2814': ['START'], 'b2815': ['START'],
                 'b0216': ['GAC', 'GAU'], 'b3797': ['CAC', 'CAU'],
                 'b4008': ['GAA', 'GAG'], 'b3798': ['CUG'],
                 'b3799': ['CCG', 'CCU', 'CCA'], 'b2652': ['AUA'],
                 'b1910': ['UGC', 'UGU'], 'b1911': ['GGC', 'GGU'],
                 'b0536': ['AGA'], 'b3658': ['UGA'], 'b1665': ['GUC', 'GUU'],
                 'b1666': ['GUC', 'GUU'], 'b3796': ['CGG'],
                 'b0744': ['GUA', 'GUG', 'GUU'], 'b0745': ['AAA', 'AAG'],
                 'b0746': ['GUA', 'GUG', 'GUU'], 'b0747': ['AAA', 'AAG'],
                 'b0743': ['AAA', 'AAG'], 'b0748': ['AAA', 'AAG'],
                 'b0749': ['AAA', 'AAG']}