from Bio import SeqIO
import pandas
from six import iteritems

from minime import util
from minime.util import dogma
from minime import *
from cobra.core import Reaction
from ecolime import ecoli_k12
from ecolime.ecoli_k12 import *
from ecolime import ribosome
import cobra
import itertools

print "Building a full ME model. Run prepare_for_minime() to make minimal " \
      "model"
macromolecules = True


def prepare_for_minime():
    global macromolecules
    macromolecules = False
    print "Creating minimal me model"


def add_transcription_reaction(me_model, TU_name, locus_ids, sequence,
                               update=True):
    """add a transcription reaction"""

    transcription = TranscriptionReaction("transcription_" + TU_name)
    transcription.transcription_data = TranscriptionData(TU_name, me_model)
    transcription.transcription_data.nucleotide_sequence = sequence
    transcription.transcription_data.RNA_products = {"RNA_" + i
                                                     for i in locus_ids}
    if not macromolecules:
        transcription.transcription_data.using_RNAP = False
    me_model.add_reaction(transcription)
    if update:
        transcription.update()
    return transcription


def create_transcribed_gene(me_model, bnum, left_pos, right_pos, seq,
                            strand, RNA_type):
    """creates a TranscribedGene object and adds it to the model"""
    gene = TranscribedGene('RNA_' + bnum)
    gene.left_pos = left_pos
    gene.right_pos = right_pos
    gene.RNA_type = RNA_type
    gene.strand = strand
    gene.nucleotide_sequence = seq

    me_model.add_metabolites([gene])
    return gene


def add_translation_reaction(me_model, bnum, amino_acid_sequence=None,
                             dna_sequence=None, update=False):
    if not amino_acid_sequence and not dna_sequence:
        print 'Transltion reactions require sequences for', bnum
    translation = TranslationReaction("translation_" + bnum)
    me_model.add_reaction(translation)

    translation.translation_data = \
        TranslationData(bnum, me_model, "RNA_" + bnum,
                        "protein_" + bnum)

    translation.translation_data.nucleotide_sequence = dna_sequence
    # translation.translation_data.get_codon_count_from_DNA(dna_sequence)
    translation.translation_data.get_last_codon_from_DNA(dna_sequence)

    amino_acid_sequence = None
    if amino_acid_sequence is not None:
        translation.translation_data.amino_acid_sequence = \
            amino_acid_sequence # .replace("U", "C")  # TODO selenocystine
    #    translation.translation_data.compute_sequence_from_DNA(dna_sequence)
    if not macromolecules:
        translation.translation_data.using_ribosome = False
    if update:
        translation.update()


def add_demand_reaction(me_model, bnum):
    warn("deprecated")

    #r = Reaction('DM_RNA'_+ bnum)
    #me_model.add_reaction(r)
    #r.reaction = 'RNA_' + bnum + ' -> '


def convert_aa_codes_and_add_charging(me_model, tRNA_aa):
    # convert amino acid 3 letter codes to metabolites
    for tRNA, aa in list(tRNA_aa.items()):
        if aa == "OTHER":
            tRNA_aa.pop(tRNA)
        elif aa == "Sec":
            #Charge with precursor to selenocysteine
            #tRNA_aa[tRNA] = me_model.metabolites.get_by_id('cys__L_c')
            tRNA_aa.pop(tRNA)
        elif aa == "Gly":
            tRNA_aa[tRNA] = me_model.metabolites.get_by_id("gly_c")
        else:
            tRNA_aa[tRNA] = \
                me_model.metabolites.get_by_id(aa.lower() + "__L_c")

    # add in all the tRNA charging reactions
    for tRNA, aa in tRNA_aa.items():
        for codon in tRNA_to_codon[tRNA]:
            tRNA_data = tRNAData("tRNA_" + tRNA + "_" + codon, me_model, aa.id,
                                 "RNA_" + tRNA, codon)
            charging_reaction = tRNAChargingReaction("charging_tRNA_" + tRNA +
                                                     "_" + codon)
            charging_reaction.tRNAData = tRNA_data
            me_model.add_reaction(charging_reaction)
            charging_reaction.update()


def build_reactions_from_genbank(me_model, gb_filename, TU_frame=None,
                                 element_types={"CDS", "rRNA", "tRNA", "ncRNA"}):

    # TODO handle special RNAse without type ('b3123')
    """create transcription and translation reactions from a genbank file

    TODO allow overriding amino acid names"""
    gb_file = SeqIO.read(gb_filename, 'gb')
    full_seq = str(gb_file.seq)

    # Determine initial amount of transcripts in model for
    # itertools.islice [start] value
    original_transcription_count = len(me_model.transcription_data)

    using_TUs = TU_frame is not None
    if not using_TUs:
        # generate a new TU frame where each mRNA gets its own TU
        TU_frame = pandas.DataFrame.from_dict(
            {"TU_" + i.qualifiers["locus_tag"][0]:
                {"start": int(i.location.start),
                 "stop": int(i.location.end),
                 "strand": "+" if i.strand == 1 else "-"}
             for i in gb_file.features if
             i.type in element_types},
            orient="index")

    tRNA_aa = {}

    # Create transcription reactions for each TU and DNA sequence.
    # RNA_products will be added so no need to update now
    for TU_id in TU_frame.index:
        sequence = dogma.extract_sequence(full_seq, TU_frame.start[TU_id],
                                          TU_frame.stop[TU_id],
                                          TU_frame.strand[TU_id])
        add_transcription_reaction(me_model, TU_id, set(), sequence,
                                   update=False)

    # Associate each feature (RNA_product) with a TU and add translation
    # reactions and demands
    for feature in gb_file.features:

        # Skip if not a gene used in ME construction
        if feature.type not in element_types or 'pseudo' in feature.qualifiers:
            continue

        # Assign values for all important gene attributes
        bnum = feature.qualifiers["locus_tag"][0]
        left_pos = int(feature.location.start)
        right_pos = int(feature.location.end)
        RNA_type = 'mRNA' if feature.type == 'CDS' else feature.type
        strand = '+' if feature.strand == 1 else '-'
        # every genetic entity gets transcription
        seq = dogma.extract_sequence(full_seq, left_pos, right_pos, strand)

        # Deal with genes that require a frameshift mutation
        frameshift_dict= ribosome.frameshift_dict
        frameshift_string = frameshift_dict.get(bnum)
        if len(seq) % 3 != 0 and frameshift_string:
            print 'Applying frameshift on %s' % bnum
            # Subtract 1 from start position to account for 0 indexing
            seq = dogma.return_frameshift_sequence(full_seq, frameshift_string)
            if strand == '-':
                seq = dogma.reverse_transcribe(seq)


        # Add translation reaction for mRNA
        if RNA_type == "mRNA":
            amino_acid_sequence = dogma.get_amino_acid_sequence_from_DNA(seq)
            add_translation_reaction(me_model, bnum, amino_acid_sequence,
                                     seq)

        # tRNA_aa = {'amino_acid':'tRNA'}
        elif feature.type == "tRNA":
            tRNA_aa[bnum] = feature.qualifiers["product"][0].split("-")[1]

        gene = create_transcribed_gene(me_model, bnum, left_pos,
                                       right_pos, seq, strand, RNA_type)

        # Add in a demand reaction for each mRNA in case the TU makes
        # multiple products and one needs a sink. If the demand reaction is
        # used, it means the mRNA doesn't count towards biomass
        demand_reaction = cobra.Reaction("DM_" + gene.id)
        me_model.add_reaction(demand_reaction)
        demand_reaction.add_metabolites(
            {gene: -1, me_model._biomass: -compute_RNA_mass(seq)})

        # associate with TU's
        parent_TU = TU_frame[
            (TU_frame.start <= left_pos + 1) & (TU_frame.stop >= right_pos) & (
            TU_frame.strand == strand)].index

        if len(parent_TU) == 0:
            if bnum not in ecoli_k12.no_TU_list:
                warn('No TU found for %s %s' % (RNA_type, bnum))
            TU_id = "TU_" + bnum
            parent_TU = [TU_id]
            add_transcription_reaction(me_model, TU_id, set(), seq, update=False)

        for TU_id in parent_TU:
            me_model.transcription_data.get_by_id(TU_id).RNA_products.add("RNA_" + bnum)

    if macromolecules:
        convert_aa_codes_and_add_charging(me_model, tRNA_aa)

    # add excised portions
    for transcription_data in itertools.islice(me_model.transcription_data,
                                               original_transcription_count,
                                               None):
        if len(transcription_data.RNA_products) == 0:
            continue
        RNA_types = set(transcription_data.RNA_types)
        if RNA_types == {"mRNA"}:
            continue
        seq = transcription_data.nucleotide_sequence
        counts = {i: seq.count(i) for i in ("A", "T", "G", "C")}
        for product_id in transcription_data.RNA_products:
            gene_seq = me_model.metabolites.get_by_id(product_id).nucleotide_sequence
            for b in ("A", "T", "G", "C"):
                counts[b] -= gene_seq.count(b)
        # excised bases
        # First base being a triphosphate will be handled by the reaction
        # producing an extra ppi during transcription. But generally, we add
        # triphosphate bases when transcribing, but excise monophosphate bases.
        monophosphate_counts = {dogma.transcription_table[k].replace("tp_c",
                                                                     "mp_c"): v
                                for k, v in iteritems(counts)}
        transcription_data.excised_bases = monophosphate_counts

    # update reactions
    for r in me_model.reactions:
        if isinstance(r, (TranscriptionReaction, TranslationReaction)):
            r.update()


def add_m_model_content(me_model, m_model, complex_metabolite_ids=[]):
    for met in m_model.metabolites:
        if met.id in complex_metabolite_ids:
            new_met = Complex(met.id)
        elif met.id.startswith("RNA"):
            new_met = TranscribedGene(met.id)
        else:
            new_met = Metabolite(met.id)
        new_met.name = met.name
        new_met.formula = met.formula
        new_met.compartment = met.compartment
        new_met.charge = met.charge
        new_met.annotation = met.annotation
        new_met.notes = met.notes
        me_model.add_metabolites(new_met)

    for reaction in m_model.reactions:
        if reaction.id.startswith("EX_") or reaction.id.startswith("DM_"):
            new_reaction = cobra.Reaction(reaction.id)
            me_model.add_reaction(new_reaction)
            new_reaction.lower_bound = reaction.lower_bound
            new_reaction.upper_bound = reaction.upper_bound
            for met, stoichiometry in iteritems(reaction.metabolites):
                new_reaction.add_metabolites(
                    {me_model.metabolites.get_by_id(met.id): stoichiometry})

        else:
            reaction_data = StoichiometricData(reaction.id, me_model)
            reaction_data.lower_bound = reaction.lower_bound
            reaction_data.upper_bound = reaction.upper_bound
            reaction_data._stoichiometry = {k.id: v for k, v
                                            in iteritems(reaction.metabolites)}


def add_generic_rRNAs(me_model):

    rRNA_type_list = ['generic_16s_rRNAs', 'generic_23s_rRNAs',
                      'generic_5s_rRNAs']

    for rRNA_type in rRNA_type_list:
        for rRNA in eval(rRNA_type):
            rRNA_id = 'RNA_' + rRNA
            me_model.add_metabolites([TranscribedGene(rRNA_type)])
            new_rxn = cobra.Reaction("rRNA_" + rRNA + '_to_generic')
            me_model.add_reaction(new_rxn)
            new_rxn.reaction = rRNA_id + ' <=> ' + rRNA_type


def add_generic_RNase(me_model, generic_flag=False):

    generic_RNase_list = eval('generic_RNase_list')
    if generic_flag:
        generic_RNase_list = [R.replace('mg2', 'generic_divalent').
                              replace('zn2', 'generic_divalent')
                              for R in generic_RNase_list]

    for cplx in generic_RNase_list:
        new_rxn = cobra.Reaction('complex_' + cplx + '_to_generic')
        me_model.add_reaction(new_rxn)
        new_rxn.reaction = cplx + ' <=> generic_RNase'


def add_modification_data(me_model, mod_id, mod_stoich, mod_enzyme=None):
    mod = ModificationData(mod_id, me_model)
    mod.stoichiometry = mod_stoich
    mod.enzyme = mod_enzyme
    return mod


def add_ribosome(me_model):
    ribosome_complex = ComplexData("ribosome", me_model)
    ribosome_components = ribosome_complex.stoichiometry
    ribosome_modifications = ribosome_complex.modifications

    add_generic_rRNAs(me_model)
    mod_dict = ribosome.ribosome_modifications
    for mod_id in mod_dict:
        mod_stoich = mod_dict[mod_id]['stoich']
        mod_enzyme = mod_dict[mod_id]['enzyme']
        num_mods = mod_dict[mod_id]['num_mods']
        mod = add_modification_data(me_model, mod_id, mod_stoich, mod_enzyme)
        ribosome_modifications[mod.id] = -num_mods

    ribosome_assembly = ribosome.ribosome_stoich
    for process in ribosome_assembly:
        for protein, amount in ribosome_assembly[process]['stoich'].items():
            #try:
            #    me_model.metabolites.get_by_id(protein)
            #except KeyError:
            #    me_model.add_metabolites([Complex(protein)])
            #    print "added", protein, "to ME model"

            ribosome_components[protein] += amount
            #try:
            #    ribosome_components[protein] += amount
            #except KeyError:
            #    ribosome_components[protein] = amount

    ribosome_complex.create_complex_formation()

def add_ribosomes_old(me_model):
    ribosome_complex = ComplexData("ribosome", me_model)
    ribosome_components = ribosome_complex.stoichiometry
    ribosome_modifications = ribosome_complex.modifications

    add_generic_rRNAs(me_model)
    mod_dict = eval('Ribosome_modifications_phase1')
    for mod_id in mod_dict:
        mod_stoich = mod_dict[mod_id]['stoich']
        mod_enzyme = mod_dict[mod_id]['enzyme']
        num_mods = mod_dict[mod_id]['num_mods']
        mod = add_modification_data(me_model, mod_id, mod_stoich, mod_enzyme)
        ribosome_modifications[mod.id] = -num_mods

    ribosome_components['generic_16s_rRNAs'] = 1

    # 30S Listed as [rpsA -rpsU], sra, [rplA-rplF],
    # rplI, [rplK-rplY], [rpmA-rpmJ]

    # 50S listed as [rplA-rplF],rplJ, rplI, rplK [rplM-rplY], [rpmA-rpmJ]

    ribosome_subunit_list = ['Ribosome_30s_proteins', 'Ribosome_50s_proteins']

    for ribosome_subunit in ribosome_subunit_list:
        protein_dict = eval(ribosome_subunit)
        for protein, amount in protein_dict.items():
            try:
                me_model.add_metabolites([Complex(protein)])
            except:
                pass
            ribosome_components[protein] = amount

    ribosome_components['mg2_c'] = 60

    # 50s reactions
    ribosome_components['generic_23s_rRNAs'] = 1
    ribosome_components['generic_5s_rRNAs'] = 1
    ribosome_components['mg2_c'] += 111

    # get ribosome ready for translation
    # ribosome_50 + ribosome_30 + trigger_factor -> rib_70
    ribosome_components['Tig_mono'] = 1

    # rib_70 + If_1 + If_3 -> rib_50_trigger_factor + rib30_if1_if3
    ribosome_components['InfA_mono'] = 1
    ribosome_components['InfC_mono'] = 1

    # 1 b3168_assumedMonomer_gtp(InfB_mono) + 1 rib_30_IF1_IF3 --> 1 rib_30_ini
    ribosome_components['InfB_mono'] = 1
    ribosome_components['gtp_c'] = 1

    # rib_30_ini + rib_50_trigger_factor -> ribsome_complex
    ribosome_complex.create_complex_formation()


def add_RNA_polymerase(me_model):
    RNAP_complex = ComplexData("RNA_Polymerase", me_model)
    RNAP_components = RNAP_complex.stoichiometry

    # for component, value in eval('RNA_polymerase_core').items():
    #    RNAP_components[component] = value
    RNAP_components['hRNAP'] = 1
    RNAP_complex.create_complex_formation()


def find_genes_within_TU(me_model, TU, RNA_pos_dict):
    loci = []
    for string_pos in RNA_pos_dict:
        pos = string_pos.split(',')
        if int(pos[0]) + 1 >= int(TU.start) and int(pos[1]) <= int(TU.stop):
            gene = RNA_pos_dict[string_pos]
            try:
                me_model.metabolites.get_by_id(gene).strand
            except AttributeError:
                print gene
            if me_model.metabolites.get_by_id(gene).strand == TU.strand:
                loci.append(RNA_pos_dict[string_pos].replace('RNA_', ''))

    return set(loci)


def should_TU_be_excised(me_model, loci):
    excise = False
    for locus in loci:
        try:
            gene = me_model.metabolites.get_by_id('RNA_'+locus)
            RNA_type = gene.RNA_type
            if RNA_type == 'tRNA' or RNA_type == 'rRNA' or locus == 'b3123':
                excise = True
        except:
            pass
    return excise


def find_existing_genome_region(me_model, left_pos, right_pos, TU_strand,
                                RNA_pos_dict, has_5prime_triphosphate):
    try:
        TU_name = RNA_pos_dict[str(left_pos)+','+str(right_pos)]
    except:
        TU_name = 'excised_TU_%s_%i_%i_%s' % \
            (TU_strand.replace('-', 'M').replace('+', 'P'), left_pos,
             right_pos, has_5prime_triphosphate)
        try:
            me_model.metabolites.get_by_id(TU_name)
        except:
            # print 'Creating new transcribed gene: ', TU_name
            new_TU = TranscribedGene(TU_name)
            new_TU.left_pos = left_pos
            new_TU.right_pos = right_pos
            new_TU.strand = TU_strand
            me_model.add_metabolites([new_TU])
            r = cobra.Reaction('DM_' + TU_name)
            me_model.add_reaction(r)
            r.reaction = TU_name + ' --> '

    return TU_name


def find_and_update_gene_at_left_pos(me_model, bnum_set, TU, gene_left_pos,
                                     excised_TU_portions,
                                     excised_TU_portion_count):
    # Check to see which rRNA, tRNA this segment represents
    for bnum in bnum_set:
        gene = me_model.metabolites.get_by_id('RNA_'+bnum)
        if gene.left_pos == gene_left_pos:
            gene_right_pos = gene.right_pos
            excised_portion_name = gene.id
            excised_gene = me_model.metabolites.get_by_id(gene.id)

    # tRNA, rRNA loses 5' triphosphate if cleaved
    if gene_left_pos != TU.start and TU.strand == '+':
        excised_gene.has_5prime_triphosphate = 'False'
    if gene_right_pos != TU.stop and TU.strand == '-':
        excised_gene.has_5prime_triphosphate = 'False'

    excised_TU_portions.append(excised_portion_name)
    excised_TU_portion_count += 1
    return gene_right_pos


def process_ends_of_TU_segment(me_model, TU, RNA_pos,
                               excised_TU_portions, RNA_pos_dict,
                               excised_TU_portion_count, right_flag):

    # (+) strain has not been sliced on the left side yet, preserving
    # triphosphate group
    if right_flag:
        right_pos = TU.stop
        left_pos = RNA_pos
        has_5prime_triphosphate = 'True' if TU.strand == '-' else 'False'
    else:
        right_pos = RNA_pos
        left_pos = TU.start
        has_5prime_triphosphate = 'False' if TU.strand == '-' else 'True'

    # Add TranscribedGene for this segment if not already in model
    TU_id = find_existing_genome_region(me_model, left_pos, right_pos,
                                        TU.strand, RNA_pos_dict,
                                        has_5prime_triphosphate)

    excised_TU_portions.append(TU_id)
    excised_TU_portion_count += 1


def find_and_add_TU_pieces(me_model, TU, left_combo_list, RNA_pos_dict,
                           bnum_set, TU_pieces):

    excised_TU_portions = []
    excised_TU_portion_count = 1

    # Look at RNA portion between TU left_pos
    # and first tRNA or rRNA
    if left_combo_list[0] > TU.start:
        RNA_left = left_combo_list[0]-1
        process_ends_of_TU_segment(me_model, TU, RNA_left, excised_TU_portions,
                                   RNA_pos_dict, excised_TU_portion_count,
                                   False)

    # iterate rest of possible segments
    while len(left_combo_list) > 0:
        gene_left_pos = left_combo_list[0]
        gene_right_pos = \
            find_and_update_gene_at_left_pos(me_model, bnum_set, TU,
                                             gene_left_pos,
                                             excised_TU_portions,
                                             excised_TU_portion_count)

        left_combo_list.pop(left_combo_list.index(gene_left_pos))

        # Deal with the last (far right) segment
        if len(left_combo_list) == 0:
            if TU.stop > gene_right_pos:
                process_ends_of_TU_segment(me_model, TU, gene_right_pos,
                                           excised_TU_portions, RNA_pos_dict,
                                           excised_TU_portion_count, True)

        # Add in excised TU for next segment, if not tRNA, rRNA
        else:
            next_left = left_combo_list[0]
            TU_id = find_existing_genome_region(me_model, gene_right_pos+1,
                                                next_left-1, TU.strand,
                                                RNA_pos_dict, 'False')
            excised_TU_portion_count += 1
            excised_TU_portions.append(TU_id)

    if len(excised_TU_portions) > 0:
        TU_pieces[TU.TU_id].append(excised_TU_portions)

    return TU_pieces


def find_all_TU_combos(me_model, TU, bnum_set, TU_seq, TU_pieces,
                       RNA_pos_dict):

    # Create list of all left strand positions of
    # tRNA, rRNA, sRNA that will be excised
    all_lefts = []
    for bnum in bnum_set:
        all_lefts.append(me_model.metabolites.get_by_id('RNA_'+bnum).left_pos)

    # Create list [0:number_of_possible_excised_portions]
    number_lefts_in_combo = range(len(all_lefts))
    number_lefts_in_combo.append(len(all_lefts))

    # Iterate through all possible combinations of gene left
    # positions that equal the number, i
    for i in number_lefts_in_combo:
        for left_combo in itertools.combinations(all_lefts, i):
            left_combo_list = list(left_combo)
            left_combo_list.sort()  # ascending by default

            if TU.TU_id not in TU_pieces:
                TU_pieces[TU.TU_id] = []

            if len(left_combo_list) == 0:
                TU_pieces[TU.TU_id].append([])
                continue

            TU_pieces = find_and_add_TU_pieces(me_model, TU, left_combo_list,
                                               RNA_pos_dict, bnum_set,
                                               TU_pieces)
    return TU_pieces


def add_excision_machinery(me_model):

    excision_types = ['rRNA_containing', 'monocistronic',
                      'polycistronic_wout_rRNA']
    for excision in excision_types:
        excision_dict = {}
        r = cobra.Reaction('combine_' + excision + '_excision_machinery')
        for machine in eval(excision):
            for ion in eval('divalent_list'):
                machine = machine.replace(ion, 'generic_divalent')
            excision_dict[Metabolite(machine)] = -1
        excision_dict[Metabolite(excision + '_excision_set')] = 1
        r.add_metabolites(excision_dict)
        me_model.add_reaction(r)

        # Add modification data objects for TU excision reactions
        rRNA_mod = ModificationData(excision + '_excision', me_model)
        rRNA_mod.stoichiometry = {'h2o_c': -1, 'h_c': 1}
        rRNA_mod.enzyme = excision + '_excision_set'


def count_and_add_excisions(me_model, pieces, transcription_data):

    tRNA_count = 0
    rRNA_count = 0
    sRNA_count = 0  # not supported

    RNA_products = set()

    for piece in pieces:
        RNA_object = me_model.metabolites.get_by_id(piece)
        if RNA_object.RNA_type == 'tRNA':
            tRNA_count += 1
        elif RNA_object.RNA_type == 'rRNA':
            rRNA_count += 1
        else:
            sRNA_count += 1
        RNA_products.add(piece)

    if rRNA_count > 0:
        transcription_data.modifications[
            'rRNA_containing_excision'] = len(pieces)-1
    elif tRNA_count == 1 and rRNA_count == 0:
        transcription_data.modifications[
            'monocistronic_excision'] = len(pieces)-1
    elif tRNA_count > 1 and rRNA_count == 0:
        transcription_data.modifications[
            'polycistronic_wout_rRNA_excision'] = len(pieces)-1
    else:  # only applies to rnpB
        transcription_data.modifications[
            'monocistronic_excision'] = len(pieces)-1
    return RNA_products


def splice_TUs(me_model, TU_pieces, generic_flag=False):

    # add_generic_RNase(me_model, generic_flag=generic_flag)
    add_excision_machinery(me_model)

    for TU, combos_of_pieces in TU_pieces.iteritems():
        for i, pieces in enumerate(combos_of_pieces):

            if len(pieces) < 1:  # no fragments to splice
                continue

            transcription = TranscriptionReaction(
                'transcription_' + TU + '_slice_' + str(i))
            transcription_data = TranscriptionData(
                TU + '_slice_' + str(i), me_model)
            full_TU_data = me_model.transcription_data.get_by_id(TU)
            transcription_data.nucleotide_sequence = \
                full_TU_data.nucleotide_sequence

            RNA_products = \
                count_and_add_excisions(me_model, pieces,
                                        transcription_data)

            transcription_data.RNA_products = RNA_products
            transcription.transcription_data = transcription_data

            me_model.add_reaction(transcription)
            transcription.update()


def add_TUs_and_translation(me_model, filename, TU_frame=None,
                            generic_flag=False):
    warn("deprecated - use build_reactions_from_genbank")

    gb_file = SeqIO.read(filename, 'gb')
    full_seq = str(gb_file.seq)
    RNA_pos_dict = build_reactions_from_genbank(me_model, gb_file,
                                                TU_frame=TU_frame)

    TU_pieces = {}
    all_rna_types = set()
    for index, TU in TU_frame.iterrows():
        seq = full_seq[TU.start:TU.stop]
        if TU.strand == '-':
            seq = util.dogma.reverse_transcribe(seq)

        loci = find_genes_within_TU(me_model, TU, RNA_pos_dict)

        excise = should_TU_be_excised(me_model, loci)

        RNA_types = {me_model.metabolites.get_by_id("RNA_" + i).RNA_type for i in loci}
        all_rna_types.add(tuple(sorted(RNA_types)))

        if not excise:
            add_transcription_reaction(me_model, TU.TU_id, loci, seq)

        else:
            TU_data = TranscriptionData(TU.TU_id, me_model)
            TU_data.nucleotide_sequence = seq
            TU_data.RNA_products = loci
            # Must splice TU to form tRNA, rRNA, sRNA (not supported)
            TU_pieces = find_all_TU_combos(me_model, TU, loci, seq, TU_pieces,
                                           RNA_pos_dict)
    print all_rna_types
    splice_TUs(me_model, TU_pieces, generic_flag=generic_flag)


def add_dummy_reactions(me_model, dna_seq, update=True):
    dummy = StoichiometricData("dummy_reaction", me_model)
    dummy.lower_bound = 0
    dummy.upper_bound = 1000
    dummy._stoichiometry = {}

    add_transcription_reaction(me_model, "dummy", {"dummy"}, dna_seq)

    me_model.add_metabolites(TranslatedGene("protein_" + "dummy"))
    add_translation_reaction(me_model, "dummy", dna_sequence=dna_seq,
                             update=update)

    complex_data = ComplexData("CPLX_dummy", me_model)
    complex_data.stoichiometry = {}
    complex_data.stoichiometry["protein_" + "dummy"] = 1
    complex_data.create_complex_formation()


def add_complex_stoichiometry_data(me_model, ME_complex_dict):
    for cplx, stoichiometry in ME_complex_dict.iteritems():
        complex_data = ComplexData(cplx, me_model)

        # stoichiometry is a defaultdict so much build as follows
        for complex, value in stoichiometry.items():
            complex_data.stoichiometry[complex] = value


def add_modication_data(me_model, mod_id, mod_dict, enzyme=None, keff=65):
    try:
        mod = me_model.modification_data.get_by_id(mod_id)
    except:
        mod = ModificationData(mod_id, me_model)
        mod.stoichiometry = mod_dict

    if not enzyme:
        mod.enzyme = enzyme
        mod.keff = keff


def add_complex_modification_data(me_model, modification_dict):
    for mod_complex_id, mod_complex_info in iteritems(modification_dict):

        unmod_complex_id, mods = mod_complex_info
        unmod_complex = me_model.complex_data.get_by_id(unmod_complex_id)

        cplx = ComplexData(mod_complex_id, me_model)
        cplx.stoichiometry = unmod_complex.stoichiometry
        cplx.translocation = unmod_complex.translocation
        cplx.chaperones = unmod_complex.chaperones

        for mod_comp, mod_count in iteritems(mods):
            mod_id = "mod_" + mod_comp
            cplx.modifications[mod_id] = -mod_count
            add_modication_data(me_model, mod_id, {mod_comp: -1})


def add_complex(me_model, ME_complex_dict,  modification_dict):

    add_complex_stoichiometry_data(me_model, ME_complex_dict)

    add_complex_modification_data(me_model, modification_dict)


def find_associated_complexes(rxnToModCplxDict, reaction_data,
                              rxn_info):
    try:
        complexes_list = rxnToModCplxDict[reaction_data.id]
    except KeyError:
        # These are orphans catalyzed by a dummy
        if reaction_data.id == "dummy_reaction" or \
                not rxn_info.is_spontaneous[reaction_data.id]:
            complexes_list = ["CPLX_dummy"]
        # These are truly spontaneous
        else:
            complexes_list = [None]
    return complexes_list


def add_and_update_metabolic_reaction(me_model, reaction_data, complex_id,
                                      reverse_flag, keff=65,
                                      update=False, create_new=True):
    if complex_id:
        complex_data = me_model.complex_data.get_by_id(complex_id)
    else:
        complex_data = None
    if not hasattr(reaction_data, "id"):
        reaction_data = me_model.stoichiometric_data.get_by_id(reaction_data)

    direction = "_REV_" if reverse_flag is True else "_FWD_"

    r = MetabolicReaction(reaction_data.id + direction + str(complex_id))
    me_model.add_reaction(r)
    r.keff = keff
    r.stoichiometric_data = reaction_data
    r.reverse = reverse_flag
    if complex_data is not None:
        r.complex_data = complex_data
    if update:
        r.update(create_new=create_new)


def add_complexes_and_rxn_data(me_model, reaction_data, complexes_list,
                               keff=65, update=False,
                               create_new=False):
    for complex_id in complexes_list:

        if reaction_data.lower_bound < 0:
            reverse_flag = True
            add_and_update_metabolic_reaction(me_model, reaction_data,
                                              complex_id, reverse_flag,
                                              keff=keff, update=update,
                                              create_new=create_new)
        if reaction_data.upper_bound > 0:
            reverse_flag = False
            add_and_update_metabolic_reaction(me_model, reaction_data,
                                              complex_id, reverse_flag,
                                              keff=keff, update=update,
                                              create_new=create_new)


def add_metabolic_reactions(me_model, reaction_data, rxnToModCplxDict,
                            rxn_info, update=False,
                            create_new=False):

    complexes_list = find_associated_complexes(rxnToModCplxDict,
                                               reaction_data, rxn_info)
    add_complexes_and_rxn_data(me_model, reaction_data, complexes_list,
                               update=update, create_new=create_new)
