from Bio import SeqIO

from minime import util
from minime import *
from cobra.core import Reaction
from ecolime.ecoli_k12 import *
import cobra
import itertools


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


def add_transcribed_gene_w_info(me_model, bnum, left_pos, right_pos, seq,
                                strand, RNA_type):

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


def add_translation_reaction(me_model, bnum, amino_acid_sequence=None,
                             dna_sequence=None, update=False):

    translation = TranslationReaction("translation_" + bnum)
    me_model.add_reaction(translation)

    translation.translation_data = \
        TranslationData(bnum, me_model, "RNA_" + bnum,
                        "protein_" + bnum)
    if amino_acid_sequence is not None:
        translation.translation_data.amino_acid_sequence = \
            amino_acid_sequence.replace("U", "C")  # TODO selenocystine
    elif dna_sequence is not None:
        translation.translation_data.compute_sequence_from_DNA(dna_sequence)

    if update:
        translation.update()


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


def build_reactions_from_genbank(me_model, gb_file, TU_frame=None):
    """create transcription and translation reactions from a genbank file

    TODO allow overriding amino acid names"""
    using_TUs = True if TU_frame is None else True

    full_seq = str(gb_file.seq)

    tRNA_aa = {}
    genome_pos_dict = {}

    for feature in gb_file.features:

        # Skip if not a gene used in ME construction
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

        # Create transcribed gene object with all important attributes
        gene = add_transcribed_gene_w_info(me_model, bnum, left_pos,
                                           right_pos, seq, strand, RNA_type)

        # Add demands for RNA if using TUs to not force translation
        if using_TUs:
            add_demand_reaction(me_model, bnum)
        else:  # If not using TUs add transcription reaction directly
            add_transcription_reaction(me_model, RNA_type + "_" + bnum,
                                       {bnum}, gene.seq)

        # Add translation reaction for all
        if feature.type == "CDS":
            try:
                amino_acid_sequence = feature.qualifiers["translation"][0]
                add_translation_reaction(me_model, bnum, amino_acid_sequence)
            except KeyError:
                continue

        # tRNA_aa = ['amino_acid':'tRNA']
        elif feature.type == "tRNA":
            tRNA_aa[bnum] = feature.qualifiers["product"][0].split("-")[1]

        genome_pos_dict[str(left_pos) + ',' + str(right_pos)] = 'RNA_' + bnum

    convert_aa_codes_and_add_charging(me_model, tRNA_aa)

    # update reactions
    for r in me_model.reactions:
        if isinstance(r, TranslationReaction):
            r.update()

    return genome_pos_dict


def add_transcription_translation_reactions(me_model, gb_filename,
                                            TU_filename=None):
    pass


def fix_id(id_str):
    return id_str.replace("_DASH_", "__")


def add_met_info(me_model, met_info, compartment_lookup, generic_ions=False):

    for met_id in met_info.index:
        fixed_id = fix_id(met_id)
        for compartment in met_info.compartment[met_id].split("AND"):
            compartment = compartment.strip()
            if compartment == "No_Compartment":
                print "Assigned %s to e" % met_id
                compartment = me_model.compartments["e"]
            new_met = Metabolite(
                fixed_id + "_" + compartment_lookup[compartment])
            new_met.name = met_info.name[met_id]
            new_met.formula = met_info.formula[met_id]
            me_model.add_metabolites(new_met)

    if generic_ions:
        compartment = "_c"
        ion_type_list = ['divalent', 'monovalent']
        for type in ion_type_list:
            new_met = cobra.Metabolite('generic_%s%s' % (type, compartment))
            new_met.name = 'Generic ' + type + ' ion'
            me_model.add_metabolites([new_met])

            for ion in eval(type + '_list'):  # load these lists from ecolime
                new_rxn = cobra.Reaction(ion + compartment + '_to_generic')
                ion_dict = {}
                new_met2 = me_model.metabolites.get_by_id(ion + compartment)
                ion_dict[new_met2] = -1
                ion_dict[new_met] = 1
                new_rxn.add_metabolites(ion_dict)
                me_model.add_reaction(new_rxn)


def add_rxn_info(me_model, rxn_info, rxn_dict):
    for rxn_id in rxn_info.index:
        reaction = StoichiometricData(rxn_id, me_model)
        reaction._stoichiometry = {fix_id(k): v
                                   for k, v in rxn_dict[rxn_id].items()}
        reaction.lower_bound = \
            -1000. if rxn_info.is_reversible[rxn_id] else 0.
        reaction.upper_bound = 1000.


def add_sources_and_sinks(me_model, sources_sinks, source_amounts,
                          compartment_lookup):

    sources_sinks.index = [fix_id(i) for i in sources_sinks.index]
    source_amounts.index = [fix_id(i) for i in source_amounts.index]

    for met in sources_sinks.index:
        met_id = met + "_" + compartment_lookup[sources_sinks.compartment[met]]
        # EX_ or DM_ + met_id
        reaction_id = sources_sinks.rxn_id[met][:3] + met_id
        reaction = cobra.Reaction(reaction_id)
        me_model.add_reaction(reaction)
        reaction.add_metabolites({me_model.metabolites.get_by_id(met_id): -1})
        # set bounds on exchanges
        if reaction.id.startswith("EX_") and met in source_amounts.index:
            reaction.lower_bound = -source_amounts.amount[met]


def add_ecoli_M_model_content(me_model, m_model_id, met_info, rxn_info,
                              rxn_dict, sources_sinks_info, source_amounts,
                              generic_ions=False):

    me_model.compartments = {"p": "Periplasm", "e": "Extra-organism",
                             "c": "Cytosol"}
    compartment_lookup = {v: k for k, v in me_model.compartments.items()}

    add_met_info(me_model, met_info, compartment_lookup,
                 generic_ions=generic_ions)

    add_rxn_info(me_model, rxn_info, rxn_dict)

    add_sources_and_sinks(me_model, sources_sinks_info, source_amounts,
                          compartment_lookup)


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


def add_ribosomes(me_model):
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

    add_generic_RNase(me_model, generic_flag=generic_flag)
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

    gb_file = SeqIO.read(filename, 'gb')
    full_seq = str(gb_file.seq)
    RNA_pos_dict = build_reactions_from_genbank(me_model, gb_file,
                                                TU_frame=TU_frame)

    TU_pieces = {}
    for index, TU in TU_frame.iterrows():
        seq = full_seq[TU.start:TU.stop]
        if TU.strand == '-':
            seq = util.dogma.reverse_transcribe(seq)

        loci = find_genes_within_TU(me_model, TU, RNA_pos_dict)

        excise = should_TU_be_excised(me_model, loci)

        if not excise:
            add_transcription_reaction(me_model, TU.TU_id, loci, seq)

        else:
            TU_data = TranscriptionData(TU.TU_id, me_model)
            TU_data.nucleotide_sequence = seq
            TU_data.RNA_products = loci
            # Must splice TU to form tRNA, rRNA, sRNA (not supported)
            TU_pieces = find_all_TU_combos(me_model, TU, loci, seq, TU_pieces,
                                           RNA_pos_dict)
    splice_TUs(me_model, TU_pieces, generic_flag=generic_flag)


def add_dummy_reactions(me_model, dna_seq):
    dummy = StoichiometricData("dummy_reaction", me_model)
    dummy.lower_bound = 0
    dummy.upper_bound = 1000
    dummy._stoichiometry = {}

    add_transcription_reaction(me_model, "dummy", {"dummy"}, dna_seq)

    me_model.add_metabolites(TranslatedGene("protein_" + "dummy"))
    add_translation_reaction(me_model, "dummy", dna_sequence=dna_seq,
                             update=True)

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


def add_complex(me_model, modifcation_dict, ME_complex_dict):

    add_complex_stoichiometry_data(me_model, ME_complex_dict)

    add_complex_modification_data(me_model, modifcation_dict)


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
