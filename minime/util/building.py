from Bio import SeqIO

from minime import util
from minime import *
from cobra.core import Reaction


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


def build_reactions_from_genbank(me_model, filename, using_TUs=False):
    """create transcription and translation reactions from a genbank file

    TOOD allow overriding amino acid names"""
    gb_file = SeqIO.read(filename, 'gb')
    full_seq = str(gb_file.seq)
    tRNA_aa = {}
    genome_pos_dict = {}
    for feature in gb_file.features:
        if feature.type == "CDS":
            bnum = feature.qualifiers["locus_tag"][0]
            seq = full_seq[feature.location.start:feature.location.end]
            gene = TranscribedGene('RNA_'+bnum)
            gene.left_pos = feature.location.start
            gene.right_pos = feature.location.end
            gene.RNA_type = 'mRNA'
            gene.has_5prime_triphosphate = 'True'
            genome_pos_dict[str(gene.left_pos) + ',' + str(gene.right_pos)] = gene.id
            me_model.add_metabolites([gene])
            if feature.strand == -1:
                gene.strand = '-'
                seq = util.dogma.reverse_transcribe(seq)
            else:
                gene.strand = '+'
            if not using_TUs:
                add_transcription_reaction(me_model, "mRNA_" + bnum, {bnum}, seq)
            try:
                amino_acid_sequence = feature.qualifiers["translation"][0]
            except KeyError:
                continue
                # translaion.translation_data.compute_sequence_from_DNA(dna_sequence)
            translation = TranslationReaction("translation_" + bnum)
            me_model.add_reaction(translation)
            r = Reaction('DM_' + bnum)
            me_model.add_reaction(r)
            r.reaction = 'RNA_' + bnum + ' -> '
            translation.translation_data = \
                TranslationData(bnum, me_model, "RNA_" + bnum,
                                "protein_" + bnum)
            translation.translation_data.amino_acid_sequence = \
                amino_acid_sequence.replace("U", "C")  # TODO selenocystine

        elif feature.type == "rRNA" or feature.type == "ncRNA":
            bnum = feature.qualifiers["locus_tag"][0]
            seq = full_seq[feature.location.start:feature.location.end]
            gene = TranscribedGene('RNA_'+bnum)
            gene.left_pos = feature.location.start
            gene.right_pos = feature.location.end
            if feature.strand == -1:
                gene.strand = '-'
            else:
                gene.strand = '+'
            gene.RNA_type = feature.type
            gene.has_5prime_triphosphate = 'True'
            genome_pos_dict[str(gene.left_pos) + ',' + str(gene.right_pos)] = gene.id
            me_model.add_metabolites([gene])
            r = Reaction('DM_' + bnum)
            me_model.add_reaction(r)
            r.reaction = 'RNA_' + bnum + ' -> '
            if feature.strand == -1:
                seq = util.dogma.reverse_transcribe(seq)
            if not using_TUs:
                add_transcription_reaction(me_model, feature.type + "_" + bnum,
                                           {bnum}, seq)

        elif feature.type == "tRNA":
            tRNA = feature
            bnum = feature.qualifiers["locus_tag"][0]
            seq = full_seq[feature.location.start:feature.location.end]
            gene = TranscribedGene('RNA_'+bnum)
            gene.left_pos = feature.location.start
            gene.right_pos = feature.location.end
            if feature.strand == -1:
                gene.strand = '-'
            else:
                gene.strand = '+'
            gene.RNA_type = feature.type
            gene.has_5prime_triphosphate = 'True'
            genome_pos_dict[str(gene.left_pos) + ',' + str(gene.right_pos)] = gene.id
            me_model.add_metabolites([gene])
            r = Reaction('DM_' + bnum)
            me_model.add_reaction(r)
            r.reaction = 'RNA_' + bnum + ' -> '
            if feature.strand == -1:
                seq = util.dogma.reverse_transcribe(seq)
            if not using_TUs:
                add_transcription_reaction(me_model, "tRNA_" + bnum, {bnum}, seq)
            tRNA_aa[bnum] = feature.qualifiers["product"][0].split("-")[1]

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

    # update reactions
    for r in me_model.reactions:
        if isinstance(r, TranslationReaction):
            r.update()
    return genome_pos_dict