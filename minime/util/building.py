from Bio import SeqIO

from minime import util
from minime import *
from cobra.core import Reaction
from ecolime.flat_files import *


def add_transcription_reaction(me_model, TU_name, locus_ids, sequence,
                               update=True):
    """add a transcription reaction"""

    transcription = TranscriptionReaction("transcription_" + TU_name)
    transcription.transcription_data = TranscriptionData(TU_name, me_model)
    transcription.transcription_data.nucleotide_sequence = sequence
    transcription.transcription_data.RNA_products = {"RNA_" + i
                                                     for i in locus_ids}
    me_model.add_reaction(transcription)
    if update:
        transcription.update()


def add_transcribed_gene(me_model, bnum, left_pos, right_pos, seq, strand, RNA_type):

    gene = TranscribedGene('RNA_'+bnum)
    gene.left_pos = left_pos
    gene.right_pos = right_pos
    gene.RNA_type = RNA_type
    gene.strand = strand

    if RNA_type != 'mRNA':
        gene.has_5prime_triphosphate = 'True'
        gene.seq = seq

    if strand == '-':
        gene.seq = util.dogma.reverse_transcribe(seq)

    me_model.add_metabolites([gene])
    return gene


def add_translation_reaction(me_model, bnum, amino_acid_sequence):

    # try:
    # except KeyError:
    #    continue
    # translaion.translation_data.compute_sequence_from_DNA(dna_sequence)
    translation = TranslationReaction("translation_" + bnum)
    me_model.add_reaction(translation)

    translation.translation_data = \
        TranslationData(bnum, me_model, "RNA_" + bnum,
                        "protein_" + bnum)
    translation.translation_data.amino_acid_sequence = \
        amino_acid_sequence.replace("U", "C")  # TODO selenocystine


def add_demand_reaction(me_model, bnum):

    r = Reaction('DM_' + bnum)
    me_model.add_reaction(r)
    r.reaction = 'RNA_' + bnum + ' -> '


def convert_aa_codes_and_add_charging(me_model, tRNA_aa):
    # convert amino acid 3 letter codes to metabolites
    for tRNA, aa in list(tRNA_aa.items()):
        if aa == "OTHER":
            tRNA_aa.pop(tRNA)
        elif aa == "Sec":
            print "TODO deal with selenocystine"
            tRNA_aa.pop(tRNA)
        elif aa == "Gly":
            tRNA_aa[tRNA] = me_model.metabolites.get_by_id("gly_c")
        else:
            tRNA_aa[tRNA] = \
                me_model.metabolites.get_by_id(aa.lower() + "__L_c")

    # add in all the tRNA charging reactions
    for tRNA, aa in tRNA_aa.items():
        tRNA_data = tRNAData("tRNA_" + tRNA, me_model, aa.id, "RNA_" + tRNA)
        charging_reaction = tRNAChargingReaction("charging_tRNA_" + tRNA)
        charging_reaction.tRNAData = tRNA_data
        me_model.add_reaction(charging_reaction)
        charging_reaction.update()


def build_reactions_from_genbank(me_model, filename, using_TUs=False):
    """create transcription and translation reactions from a genbank file

    TODO allow overriding amino acid names"""

    gb_file = SeqIO.read(filename, 'gb')
    full_seq = str(gb_file.seq)
    tRNA_aa = {}
    genome_pos_dict = {}

    for feature in gb_file.features:

        add_list = ['CDS', 'rRNA', 'tRNA', 'ncRNA']
        if feature.type not in add_list:
            continue

        # Assign values for all important gene attributes
        bnum = feature.qualifiers["locus_tag"][0]
        left_pos = feature.location.start
        right_pos = feature.location.end
        seq = full_seq[feature.location.start:feature.location.end]
        RNA_type = 'mRNA' if feature.type == 'CDS' else feature.type
        strand = '+' if feature.strand == 1 else '-'

        gene = add_transcribed_gene(me_model, bnum, left_pos, right_pos, seq, strand, RNA_type)

        if using_TUs:
            add_demand_reaction(me_model, bnum)
        else:
            add_transcription_reaction(me_model, RNA_type + "_" + bnum, {bnum}, gene.seq)

        if feature.type == "CDS":
            try:
                amino_acid_sequence = feature.qualifiers["translation"][0]
                add_translation_reaction(me_model, bnum, amino_acid_sequence)
            except KeyError:
                continue

        elif feature.type == "tRNA":
            tRNA_aa[bnum] = feature.qualifiers["product"][0].split("-")[1]

        genome_pos_dict[str(left_pos) + ',' + str(right_pos)] = 'RNA_' + bnum

    convert_aa_codes_and_add_charging(me_model, tRNA_aa)

    # update reactions
    for r in me_model.reactions:
        if isinstance(r, TranslationReaction):
            r.update()

    return genome_pos_dict