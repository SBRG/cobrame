from __future__ import division, absolute_import, print_function

from six import iteritems

from .dogma import transcription_table

# Mass of amino acids (Daltons) with H and OH removed from either end.
# This is used to compute protein masses because a water is removed
# during formation of the amide bond.
amino_acid_no_h2o = {
    'ala__L_c': 71.079,
    'arg__L_c': 157.197,
    'asn__L_c': 114.104,
    'asp__L_c': 114.08,
    'cys__L_c': 103.145,
    'gln__L_c': 128.131,
    'glu__L_c': 128.107,
    'gly_c': 57.052,
    'his__L_c': 137.142,
    'ile__L_c': 113.16,
    'leu__L_c': 113.16,
    'lys__L_c': 129.183,
    'met__L_c': 131.199,
    'phe__L_c': 147.177,
    'pro__L_c': 97.117,
    'ser__L_c': 87.078,
    'thr__L_c': 101.105,
    'trp__L_c': 186.214,
    'tyr__L_c': 163.176,
    'val__L_c': 99.133}

# Mass of the RNA nucleotides with the diphosphate removed.
rna_no_ppi = {
    'gtp_c': 344.2,
    'utp_c': 305.159,
    'ctp_c': 304.175,
    'atp_c': 328.201}

# Mass of the DNA nucleotides with the diphosphate removed.
dna_mw_no_ppi = {'datp': 312.202,
                 'dctp': 286.16,
                 'dgtp': 328.201,
                 'dttp': 303.187}


def compute_protein_mass(amino_acid_count):
    """compute protein mass in kDa from amino acid count

    amino_acid_count: {amino_acid: number}

    """
    protein_mass = sum(amino_acid_no_h2o[aa] * count
                       for aa, count in iteritems(amino_acid_count))
    protein_mass += 18.015  # one water not removed
    return protein_mass / 1000


def compute_rna_mass(dna_sequence, excised_bases=None):
    """compute RNA mass in kDA from nucleotide count

    nucleotide_count: {nucleotide: number}

    """
    if not excised_bases:
        excised_bases = {}

    nuc_count = {transcription_table[i]: dna_sequence.count(i)
                 for i in set(dna_sequence)}
    nuc_count["gtp_c"] -= excised_bases.get("gmp_c", 0)
    nuc_count["utp_c"] -= excised_bases.get("ump_c", 0)
    nuc_count["ctp_c"] -= excised_bases.get("cmp_c", 0)
    nuc_count["atp_c"] -= excised_bases.get("amp_c", 0)
    rna_mass = sum(rna_no_ppi[nuc] * count
                   for nuc, count in iteritems(nuc_count))
    if sum(excised_bases.values()) > 0:
        rna_mass += 174.951262  # 5' has 3 phosphates
    return rna_mass / 1000
