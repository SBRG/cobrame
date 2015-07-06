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


def add_transcribed_gene(me_model, bnum, left_pos, right_pos, seq, strand,
                         RNA_type):

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
        gene = add_transcribed_gene(me_model, bnum, left_pos,
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
            new_met = cobra.Metabolite('generic_%s_%s' % (type, compartment))
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
        if reaction.id.startswith("EX_") and met_id in source_amounts.index:
            reaction.lower_bound = -source_amounts.amount[met_id]


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

    rRNA_type_list = ['generic_16s_rRNAs', 'generic_23s_rRNAs', 'generic_5s_rRNAs']

    for rRNA_type in rRNA_type_list:
        for rRNA in eval(rRNA_type):
            rRNA_id = 'RNA_' + rRNA
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

    ribosome_components['generic_16'] = 1

    # 30S Listed as [rpsA -rpsU], sra, [rplA-rplF],
    # rplI, [rplK-rplY], [rpmA-rpmJ]

    # 50S listed as [rplA-rplF],rplJ, rplI, rplK [rplM-rplY], [rpmA-rpmJ]

    ribosome_subunit_list = ['Ribosome_30s_proteins', 'Ribosome_50s_proteins']

    for ribosome_subunit in ribosome_subunit_list:
        protein_dict = eval(ribosome_subunit)
        for protein, amount in protein_dict.items():
            ribosome_components[protein] = amount

    ribosome_components['mg2_c'] = 60

    # 50s reactions
    ribosome_components['generic_23s'] = 1
    ribosome_components['generic_5s'] = 1
    ribosome_components['mg2_c'] += 111

    # get ribosome ready for translation
    # ribosome_50 + ribosome_30 + trigger_factor -> rib_70
    ribosome_components['Tig_mono'] = 1

    # rib_70 + If_1 + If_3 -> rib_50_trigger_factor + rib30_if1_if3
    ribosome_components['InfA_mono'] = 1
    ribosome_components['InfC_mono'] = 1

    # 1 b3168_assumedMonomer_gtp (InfB_mono) + 1 rib_30_IF1_IF3 --> 1 rib_30_ini
    ribosome_components['InfB_mono'] = 1
    ribosome_components['gtp_c'] = 1

    # rib_30_ini + rib_50_trigger_factor -> ribsome_complex
    ribosome_complex.create_complex_formation()


def add_RNA_polymerase(me_model):
    RNAP_complex = ComplexData("RNA_Polymerase", me_model)
    RNAP_components = RNAP_complex.stoichiometry

    #for component, value in eval('RNA_polymerase_core').items():
    #    RNAP_components[component] = value
    RNAP_components['hRNAP'] = 1
    RNAP_complex.create_complex_formation()


def find_genes_within_TU(me_model, TU, RNA_pos_dict):
    loci = []
    for string_pos in RNA_pos_dict:
        pos = string_pos.split(',')
        if int(pos[0]) + 1 >= int(TU.start) and int(pos[1]) <= int(TU.stop):
            if me_model.metabolites.get_by_id(RNA_pos_dict[string_pos]).strand == TU.strand:
                loci.append(RNA_pos_dict[string_pos].replace('RNA_', ''))

    return set(loci)


def should_TU_be_excised(me_model, loci):
    excise = False
    for locus in loci:
        try:
            gene = me_model.metabolites.get_by_id('RNA_'+locus)
            if gene.RNA_type == 'tRNA' or ('RNA_'+locus).RNA_type == 'rRNA' or locus == 'b3123':
                excise = True
        except:
            pass
    return excise


def find_existing_genome_region(me_model, left_pos, right_pos, TU_strand,
                                RNA_pos_dict, has_5prime_triphosphate):
    try:
        TU_name = RNA_pos_dict[str(left_pos)+','+str(right_pos)]
    except:
        TU_name = 'excised_TU_%s_%i_%i_%s' % (TU_strand.replace('-','M').replace('+','P'),
                                                  left_pos, right_pos, has_5prime_triphosphate)
        try:
            me.metabolites.get_by_id(TU_name)
        except:
            #print 'Creating new transcribed gene: ', TU_name
            new_TU = TranscribedGene(TU_name)
            new_TU.left_pos = left_pos
            new_TU.right_pos = right_pos
            new_TU.strand = TU_strand
            me_model.add_metabolites([new_TU])
            r = cobra.Reaction('DM_' + TU_name)
            me_model.add_reaction(r)
            r.reaction = TU_name + ' --> '

    return TU_name


def find_and_update_gene_at_left_pos(me_model, bnum_set, TU_left,
                                     gene_left_pos, TU_strand,
                                     excised_TU_portions,
                                     excised_TU_portion_count):
    # Check to see which rRNA, tRNA this segment represents
    for bnum in bnum_set:
        gene = me_model.metabolites.get_by_id('RNA_'+bnum)
        if gene.left_pos == gene_left_pos:
            gene_right_pos = gene.right_pos
            excised_portion_name = gene_id
            excised_gene = me_model.metabolites.get_by_id(gene_id)

    # tRNA, rRNA loses 5' triphosphate if cleaved
    if gene_left_pos != TU_left and TU_strand == '+':
        excised_gene.has_5prime_triphosphate = 'False'
    if gene_right_pos != TU_right and TU_strand == '-':
        excised_gene.has_5prime_triphosphate = 'False'

    excised_TU_portions.append(excised_portion_name)
    excised_TU_portion_count += 1
    return gene_right_pos


def process_ends_of_TU_segment(me_model, TU_left, RNA_left_pos, TU_strand,
                                excised_TU_portions,
                                excised_TU_portion_count, right_flag):

    # (+) strain has not been sliced on the left side yet, preserving
    # triphosphate group
    if right_flag:
        has_5prime_triphosphate = 'True' if TU_strand == '-' else 'False'
    else:
        has_5prime_triphosphate = 'False' if TU_strand == '-' else 'True'


    # Add TranscribedGene for this segment if not already in model
    TU_id = find_existing_genome_region(me_model, TU_left, RNA_left_pos,
                                        TU_strand, has_5prime_triphosphate)

    excised_TU_portions.append(TU_id)
    excised_TU_portion_count += 1


def transcribe_all_TU_combos(me_model, TU, TU_left, TU_right, TU_strand,
                             bnum_set, TU_seq, TU_pieces, RNA_pos_dict):

    # Create list of all left strand positions of tRNA, rRNA, sRNA that will be excised
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

            if TU not in TU_pieces:
                TU_pieces[TU] = []

            if len(left_combo_list) == 0:
                TU_pieces[TU].append([])
                continue

            excised_TU_portions = []
            excised_TU_portion_count = 1

            # Look at RNA portion between TU left_pos
            # and first tRNA or rRNA
            if left_combo_list[0] > TU_left:
                process_ends_of_TU_segment(me_model, TU_left, left_combo_list[0]-1,
                                            TU_strand, excised_TU_portions,
                                            excised_TU_portion_count, False)

            # iterate rest of possible segments
            while len(left_combo_list) > 0:
                gene_left_pos = left_combo_list[0]
                gene_right_pos = \
                    find_and_update_gene_at_left_pos(me_model, bnum_set, TU_left,
                                                     gene_left_pos, TU_strand,
                                                     excised_TU_portions,
                                                     excised_TU_portion_count)

                left_combo_list.pop(left_combo_list.index(gene_left_pos))

                # Deal with the last (far right) segment
                if len(left_combo_list) == 0:
                    if TU_right > gene_right_pos:

                        if left_combo_list[0] > TU_left:
                                process_ends_of_TU_segment(me_model, TU_left, left_combo_list[0]-1,
                                            TU_strand, excised_TU_portions,
                                            excised_TU_portion_count, True)
                        # (-) strain has not been sliced on the right side yet, preserving
                        # triphosphate group
                        if TU_strand == '-':
                            has_5prime_triphosphate = 'True'
                        else:
                            has_5prime_triphosphate = 'False'

                        TU_id = find_existing_genome_region(me_model, gene_right_pos+1, TU_right,
                                                            TU_strand, RNA_pos_dict,
                                                            has_5prime_triphosphate)

                        excised_TU_portions.append(TU_id)
                        excised_TU_portion_count += 1

                # Add in excised TU for next segment, if not tRNA, rRNA
                else:
                    next_left = left_combo_list[0]
                    TU_id = find_existing_genome_region(me_model, gene_right_pos+1,
                                                        next_left-1, TU_strand, 'False')
                    excised_TU_portion_count += 1
                    excised_TU_portions.append(TU_id)

            if len(excised_TU_portions) > 0:
                TU_pieces[TU].append(excised_TU_portions)


def add_TUs_and_translation(me_model, gb_file, TU_frame, RNA_pos_dict):
    full_seq = str(gb_file.seq)
    TU_pieces = {}

    for index, TU in TU_frame.iterrows():
        seq = full_seq[TU.start:TU.stop]
        if TU.strand == '-':
            seq = util.dogma.reverse_transcribe(seq)

        loci = find_genes_within_TU(me_model, TU, RNA_pos_dict)

        excise = should_TU_be_excised(me_model, loci)

        if not excise:
            add_transcription_reaction(me_model, index, loci, seq)

        else:
            TU_data = TranscriptionData(index, me_model)
            TU_data.nucleotide_sequence = seq
            TU_data.RNA_products = loci
            # Must splice TU to form tRNA, rRNA, sRNA (not supported)
            transcribe_all_TU_combos(me_model, index, TU.start, TU.stop,
                                     TU.strand, loci, seq, TU_pieces,
                                     RNA_pos_dict)