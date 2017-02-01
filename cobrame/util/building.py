from Bio import SeqIO
import pandas
from six import iteritems

from cobrame import util
from cobrame.util import dogma
from cobrame import *
from cobra.core import Reaction
from ecolime.ecoli_k12 import *
import cobra
import itertools


def add_transcription_reaction(me_model, TU_name, locus_ids, sequence,
                               update=True):
    """
    Create TranscriptionReaction object and add it to ME-Model.
    This includes the necessary transcription data.

    Args:
        me_model: cobra.model.MEModel
            The MEModel object to which the reaction will be added

        TU_name: String
            ID of TU being transcribed.
            The TranscriptionReaction will be added as "transcription_+TU_name"
            The TranscriptionData will be added as just "TU_name"

        locus_ids: Set
            Set of locus IDs that the TU transcribes

        sequence: String
            Nucleotide sequence of the TU.

        update: Boolean
            If True, use TranscriptionReaction's update function to update and
            add reaction stoichiometry

    Returns: cobra.core.Reaction.TranscriptionReaction
        Return TranscriptionReaction for the TU
    """

    transcription = TranscriptionReaction("transcription_" + TU_name)
    transcription.transcription_data = TranscriptionData(TU_name, me_model)
    transcription.transcription_data.nucleotide_sequence = sequence
    transcription.transcription_data.RNA_products = {"RNA_" + i
                                                     for i in locus_ids}

    me_model.add_reaction(transcription)
    if update:
        transcription.update()
    return transcription


def create_transcribed_gene(me_model, locus_id, left_pos, right_pos, seq,
                            strand, RNA_type):
    """
     Creates a TranscribedGene metabolite object and adds it to the ME-model

    Args:
        me_model: cobra.model.MEModel
            The MEModel object to which the reaction will be added

        locus_id: String
            Locus ID of RNA product.
            The TranscribedGene will be added as "RNA_ + locus_id"

        left_pos: Integer
            Left position of gene on the sequence of the (+) strain

        right_pos: Integer
            Right position of gene on the sequence of the (+) strain

        seq: String
            Nucleotide sequence of RNA product.
            Amino acid sequence, codon counts, etc. will be calculated based on
            this string.

        strand: String
            (+) if the RNA product is on the leading strand
            (-) if the RNA product is on the complementary strand

        RNA_type: String
            Type of RNA of the product.
            tRNA, rRNA, or mRNA
            Used for determining how RNA product will be processed.

    Returns:
        cobra.core.Metabolite.TranscribedGene object for the RNA product
    """
    gene = TranscribedGene('RNA_' + locus_id)
    gene.left_pos = left_pos
    gene.right_pos = right_pos
    gene.RNA_type = RNA_type
    gene.strand = strand
    gene.nucleotide_sequence = seq

    me_model.add_metabolites([gene])
    return gene


def add_translation_reaction(me_model, locus_id, dna_sequence=None,
                             update=False):
    """
    Creates and adds a TranslationReaction to the ME-model as well as the
    associated TranslationData

    A dna_sequence  is required in order to add a TranslationReaction to the
    ME-model

    Args:
        me_model: cobra.model.MEModel
            The MEModel object to which the reaction will be added

        locus_id: String
            Locus ID of RNA product.
            The TranslationReaction will be added as "translation_ + locus_id"
            The TranslationData will be added as "locus_id"

        dna_sequence: String
            DNA sequence of the RNA product. This string should be reverse
            transcribed if it originates on the complement strand.

        update: Boolean
            If True, use TranslationReaction's update function to update and add
            reaction stoichiometry

        terminator_dict: Dict
            {stop_codon: release_factor}

            Used to determine which ProcessData.SubReaction to add to the
            TranslationReaction to account for termination of peptide
    """

    # Create TranslationData
    translation_data = TranslationData(locus_id, me_model, "RNA_" + locus_id,
                                       "protein_" + locus_id)
    translation_data.nucleotide_sequence = dna_sequence

    # Add RNA to model if it doesn't exist
    if "RNA_" + locus_id not in me_model.metabolites:
        RNA = TranscribedGene('RNA_c')
        RNA.nucleotide_sequence = dna_sequence
        me_model.add_metabolites(RNA)

    # Create and add TranslationReaction with TranslationData
    translation_reaction = TranslationReaction("translation_" + locus_id)
    me_model.add_reaction(translation_reaction)
    translation_reaction.translation_data = translation_data

    if update:
        translation_reaction.update()


def convert_aa_codes_and_add_charging(me_model, tRNA_aa, tRNA_to_codon,
                                      verbose=True):
    """
    Adds tRNA charging reactions for all tRNAs in ME-model

    Args:
        me_model: cobra.model.MEModel
            The MEModel object to which the reaction will be added

        tRNA_aa: Dict
            Dictionary of tRNA locus ID to 3 letter codes of the amino acid
            that the tRNA contributes

            {locus_id: amino_acid_3_letter_code}

        verbose: Boolean
            If True, display metabolites that were not previously added to the
            model and were thus added when creating charging reactions
    """
    # convert amino acid 3 letter codes to metabolites
    for tRNA, aa in list(tRNA_aa.items()):
        if aa == "OTHER":
            tRNA_aa.pop(tRNA)
        elif aa == "Sec":
            # Charge with precursor to selenocysteine
            tRNA_aa[tRNA] = me_model.metabolites.get_by_id('cys__L_c')
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
            charging_reaction.update(verbose=verbose)


def build_reactions_from_genbank(me_model, gb_filename, TU_frame=None,
                                 element_types={'CDS', 'rRNA', 'tRNA', 'ncRNA'},
                                 verbose=True, frameshift_dict={},
                                 tRNA_to_codon={}, update=True):

    # TODO handle special RNAse without type ('b3123')
    """Creates and adds transcription and translation reactions using genomic
     information from the organism's genbank file. Adds in the basic
     requirements for these reactions. Organism specific components are added
     ...

    Args:
        me_model: cobra.model.MEModel
            The MEModel object to which the reaction will be added

        gb_filename: String
            Local name of the genbank file that will be used for ME-model
            construction

        TU_frame: pandas.DataFrame
            DataFrame with indexes of the transcription unit name and columns
            containing the transcription unit starting and stopping location on
            the genome and whether the transcription unit is found on the
            main (+) strand or complementary (-) strand.

            If no transcription unit DataFrame is passed into the function,
            transcription units are added corresponding to each transcribed
            gene in the genbank file.

        element_types: Set
            Transcription reactions will be added to the ME-model for all RNA
            feature.types in this set. This uses the nomenclature of the
            genbank file (gb_filename)


        verbose: Boolean
            If True, display metabolites that were not previously added to the
            model and were thus added when creating charging reactions

        translation_terminators: Dict
            {stop_codon: release_factor}

            Used to determine which ProcessData.SubReaction to add to the
            TranslationReaction to account for termination of peptide

        frameshift_dict: Dict
            {locus_id: genome_position_of_TU}

            If a locus_id is in the frameshift_dict, update it's nucleotide
            sequence to account of the frameshift

    """

    # Load genbank file and extract DNA sequence
    gb_file = SeqIO.read(gb_filename, 'gb')
    full_seq = str(gb_file.seq)

    # Dictionary of tRNA locus ID to the 3 letter code for the amino acid it
    # contributes
    tRNA_aa = {}

    # If no TU_frame is provided generate a new TU frame where each mRNA gets
    # its own TU
    using_TUs = TU_frame is not None
    if not using_TUs:
        TU_frame = pandas.DataFrame.from_dict(
            {"TU_" + i.qualifiers["locus_tag"][0]:
                {"start": int(i.location.start),
                 "stop": int(i.location.end),
                 "strand": "+" if i.strand == 1 else "-"}
             for i in gb_file.features if
             i.type in element_types},
            orient="index")

    # Create transcription reactions for each TU and DNA sequence.
    # RNA_products will be added so no need to update now
    for TU_id in TU_frame.index:
        # subtract 1 from TU start site to account for 0 indexing
        sequence = dogma.extract_sequence(full_seq, TU_frame.start[TU_id]-1,
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

        # ---- Assign values for all important gene attributes ----
        bnum = feature.qualifiers["locus_tag"][0]
        left_pos = int(feature.location.start)
        right_pos = int(feature.location.end)
        RNA_type = 'mRNA' if feature.type == 'CDS' else feature.type
        strand = '+' if feature.strand == 1 else '-'
        seq = dogma.extract_sequence(full_seq, left_pos, right_pos, strand)

        # ---- Add gene metabolites and apply frameshift mutations----
        frameshift_string = frameshift_dict.get(bnum)
        if len(seq) % 3 != 0 and frameshift_string:
            print('Applying frameshift on %s' % bnum)
            seq = dogma.return_frameshift_sequence(full_seq, frameshift_string)
            if strand == '-':
                seq = dogma.reverse_transcribe(seq)

        # Add TranscribedGene metabolite
        gene = create_transcribed_gene(me_model, bnum, left_pos,
                                       right_pos, seq, strand, RNA_type)

        # ---- Add translation reaction for mRNA ----
        if RNA_type == "mRNA":
            add_translation_reaction(me_model, bnum, dna_sequence=seq)

        # ---- Create dict to use for adding tRNAChargingReactions ----
        # tRNA_aa = {'amino_acid':'tRNA'}
        elif RNA_type == "tRNA":
            tRNA_aa[bnum] = feature.qualifiers["product"][0].split("-")[1]

        # ---- Add in a demand reaction for each mRNA ---
        # This is in case the TU makes multiple products and one needs a sink.
        # If the demand reaction is used, it means the mRNA doesn't count
        # towards biomass
        demand_reaction = cobra.Reaction("DM_" + gene.id)
        me_model.add_reaction(demand_reaction)
        demand_reaction.add_metabolites({gene: -1})

        # mRNA biomass is handled during translation
        if RNA_type == 'tRNA':
            demand_reaction.add_metabolites({
                me_model._tRNA_biomass: -compute_RNA_mass(seq)})
        elif RNA_type == 'rRNA':
            demand_reaction.add_metabolites({
                me_model._rRNA_biomass: -compute_RNA_mass(seq)})
        elif RNA_type == 'ncRNA':
            demand_reaction.add_metabolites({
                me_model._ncRNA_biomass: -compute_RNA_mass(seq)})
        elif RNA_type == 'mRNA':
            demand_reaction.add_metabolites({
                me_model._mRNA_biomass: -compute_RNA_mass(seq)})

        # ---- Associate TranscribedGene to a TU ----
        parent_TU = TU_frame[
            (TU_frame.start - 1 <= left_pos) & (TU_frame.stop >= right_pos) & (
            TU_frame.strand == strand)].index

        if len(parent_TU) == 0:
            if verbose:
                warn('No TU found for %s %s' % (RNA_type, bnum))
            TU_id = "TU_" + bnum
            parent_TU = [TU_id]
            add_transcription_reaction(me_model, TU_id, set(), seq,
                                       update=False)

        for TU_id in parent_TU:
            me_model.transcription_data.get_by_id(TU_id).RNA_products.add(
                    "RNA_" + bnum)

    convert_aa_codes_and_add_charging(me_model, tRNA_aa, tRNA_to_codon,
                                      verbose=verbose)

    if update:
        for r in me_model.reactions:
            if isinstance(r, (TranscriptionReaction, TranslationReaction)):
                r.update()


def add_m_model_content(me_model, m_model, complex_metabolite_ids=[]):
    """
    Add metabolite and reaction attributes to me_model from m_model. Also
    creates StoichiometricData objects for each reaction in m_model, and adds
    reactions directly to me_model if they are exchanges or demands.

    Args:
        me_model: cobra.model.MEModel
            The MEModel object to which the content will be added

        m_model: cobra.model
            The m_model which will act as the source of metabolic content for
            MEModel

        complex_metabolite_ids: list
            List of complexes which are 'metabolites' in the m-model reaction
            matrix, but should be treated as complexes

    """
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


def add_dummy_reactions(me_model, dna_seq, update=True):
    dummy = StoichiometricData("dummy_reaction", me_model)
    dummy.lower_bound = 0
    dummy.upper_bound = 1000
    dummy._stoichiometry = {'CPLX_dummy': -1}

    create_transcribed_gene(me_model, 'dummy', 0, len(dna_seq), dna_seq, '+',
                            'mRNA')
    add_transcription_reaction(me_model, "RNA_dummy", {"dummy"}, dna_seq)
    me_model.add_metabolites(TranslatedGene("protein_" + "dummy"))
    add_translation_reaction(me_model, "dummy", dna_sequence=dna_seq,
                             update=update)
    try:
        complex_data = ComplexData("CPLX_dummy", me_model)
    except ValueError:
        warn('CPLX_dummy already in model')
        complex_data = me_model.complex_data.get_by_id('CPLX_dummy')
    complex_data.stoichiometry = {"protein_dummy": 1}
    if update:
        complex_data.create_complex_formation()


def add_complex_stoichiometry_data(me_model, ME_complex_dict):
    warn('deprecated')
    for cplx, stoichiometry in ME_complex_dict.iteritems():
        complex_data = ComplexData(cplx, me_model)

        # stoichiometry is a defaultdict so much build as follows
        for complex, value in stoichiometry.items():
            complex_data.stoichiometry[complex] = value


def add_complex_to_model(me_model, complex_id, complex_stoichiometry,
                         complex_modifications={}):
    """
    Adds ComplexData to the model for a given complex.

    Args:
        me_model: cobrame.MEModel

        complex_id: string
            ID of the complex and thus the model ComplexData

        complex_stoichiometry: dict
            {complex_id: {protein_<locus_tag>: stoichiometry}}

        complex_modifications: dict
            {modification_id: stoichiometry}

    """

    complex_data = ComplexData(complex_id, me_model)
    # must add update stoichiometry one by one since it is a defaultdict
    for metabolite, value in complex_stoichiometry.items():
        complex_data.stoichiometry[metabolite] += value
    for modification, value in complex_modifications.items():
        complex_data.modifications[modification] = value


def add_modification_data(me_model, modification_id,
                          modification_stoichiometry,
                          modification_enzyme=None, verbose=True):
    """
    Creates a ModificationData object for each modification defined by the
    function inputs.

    It's assumed every complex modification occurs spontaneously, unless

    If a modification uses an enzyme this can be updated after the
    ModificationData object is already created


    Args:
        me_model: class:cobrame.MEModel
        metabolite: str
            ID of the metabolite that the enzyme complex is being modified by

    """

    if modification_id in me_model.modification_data:
        if verbose:
            warn('Modification (%s) already in model' % modification_id)
        else:
            pass
    else:
        modification_data = ModificationData(modification_id, me_model)
        modification_data.stoichiometry = modification_stoichiometry
        modification_data.enzyme = modification_enzyme


def add_model_complexes(me_model, complex_stoichiometry_dict,
                        complex_modification_dict, verbose=True):
    """
    Construct ComplexData for complexes into MEModel from its subunit
    stoichiometry, and a dictionary of its modification metabolites.

    It is assumed that each modification adds one equivalent of the
    modification metabolite. Multiple


    Intended to be used as a function for large-scale complex addition.

    For adding individual ComplexData objects, use add_complex_to_model

    Args:
        me_model: cobrame.MEModel

        complex_stoichiometry_dict: dict
            {unmodified_complex_id: {protein_<locus_tag>: stoichiometry}}

        complex_modification_dict: dict
            {modified_complex_id:{core_enzyme: unmodified_complex_id,
                                 'modifications: {mod_metabolite:
                                                  stoichiometry}}}

    """
    for complex_id, stoichiometry in complex_stoichiometry_dict.items():
        add_complex_to_model(me_model, complex_id, stoichiometry, {})

    for modified_complex_id, info in complex_modification_dict.items():
        modification_dict = {}
        for metabolite, number in info['modifications'].items():
            modification_id = 'mod_' + metabolite
            add_modification_data(me_model, modification_id, {metabolite: -1},
                                  verbose=verbose)
            # stoichiometry of modification determined in
            # modification_data.stoichiometry
            modification_dict[modification_id] = abs(number)

        core_enzyme = \
            complex_modification_dict[modified_complex_id]['core_enzyme']
        stoichiometry = complex_stoichiometry_dict[core_enzyme]

        add_complex_to_model(me_model, modified_complex_id, stoichiometry,
                             complex_modifications=modification_dict)


def add_metabolic_reaction_to_model(me_model, stoichiometric_data_id,
                                    directionality, complex_id=None,
                                    spontaneous=False, update=False,
                                    keff=65):
    """
    Creates and add a MetabolicReaction to a MEModel.


    Args:
        me_model:
            MEModel that the MetabolicReaction will be added to

        stoichiometric_data_id: string
            ID of the StoichiometricData for the reaction being added

        directionality: string
            Forward: Add reaction that occurs in the forward direction
            Reverse: Add reaction that occurs in the reverse direction

        complex_id: string or None
            ID of the ComplexData for the enzyme that catalyze the reaction
            being added.

        spontaneous: boolean
            If True and complex_id='' add reaction as spontaneous reaction
            If False and complex_id='' add reaction as orphan (CPLX_dummy
            catalyzed)

    """
    # Get stoichiometric data for reaction being added
    try:
        stoichiometric_data = me_model.stoichiometric_data.get_by_id(stoichiometric_data_id)
    except KeyError:
        raise Exception("Stoichiometric data for %s has not been added to"
                        " model" % stoichiometric_data_id)

    # Get complex data and id based on arguments passed into function
    if type(complex_id) == str:
        complex_data = me_model.complex_data.get_by_id(complex_id)
    elif complex_id is None and spontaneous is True:
        complex_id = "SPONT"
        complex_data = None
    elif complex_id is None and spontaneous is False:
        complex_id = "CPLX_dummy"
        try:
            complex_data = me_model.complex_data.get_by_id(complex_id)
        except KeyError:
            raise Exception("CPLX_dummy must be added to complex data to add"
                            "orphan reactions")
    else:
        raise ValueError("Complex id (%s) must be a string or None" %
                         str(complex_id))

    if directionality.lower() == 'forward':
        direction = "_FWD_"
        reverse_flag = False
    elif directionality.lower() == 'reverse':
        direction = "_REV_"
        reverse_flag = True
    else:
        raise NameError("Reaction direction must be 'forward' or 'reverse'")

    r = MetabolicReaction(stoichiometric_data_id + direction + complex_id)
    me_model.add_reaction(r)
    r.keff = keff
    r.stoichiometric_data = stoichiometric_data
    r.reverse = reverse_flag
    if complex_data is not None:
        r.complex_data = complex_data
    if update:
        r.update(verbose=True)


def add_reactions_from_stoichiometric_data(me_model, rxn_to_cplx_dict,
                                           rxn_info_frame, update=False,
                                           keff=65):
    """
    Creates and adds MetabolicReaction for all StoichiometricData in model.

    Intended for use when adding all reactions from stoichiometric data for the
    first time.

    For adding an individual reaction use add_metabolic_reaction_to_model()

    Args:
        me_model:
            MEModel that the MetabolicReaction will be added to

        rxn_to_cplx_dict: Dict
            {StoichiometricData.id: catalytic_enzyme_id}


        rxn_info_frame: pandas.Dataframe
            Contains the ids, names and reversibility for each reaction in the
            metabolic reaction matrix as well as whether the reaction is
            spontaneous
    """

    for reaction_data in me_model.stoichiometric_data:

        try:
            spontaneous_flag = rxn_info_frame.is_spontaneous[reaction_data.id]
        except KeyError:
            spontaneous_flag = 0
            warn("(%s) not in rxn_info_frame assumed nonspontaneous" %
                 reaction_data.id)

        if spontaneous_flag == 1:
            spontaneous = True
        elif spontaneous_flag == 0:
            spontaneous = False
        else:
            raise Exception("is_spontaneous must be '1' or '0'")

        # Reactions can be catalyzed by multiple isozymes so retrieve list of
        # complexes that catalyze the reaction
        try:
            complexes_list = rxn_to_cplx_dict[reaction_data.id]
        except KeyError:
            complexes_list = [None]

        # Add metabolic reactions for each isozyme
        for complex_id in complexes_list:
            directionality_list = []
            if reaction_data.lower_bound < 0:
                directionality_list.append('reverse')
            if reaction_data.upper_bound > 0:
                directionality_list.append('forward')
            for directionality in directionality_list:
                add_metabolic_reaction_to_model(me_model, reaction_data.id,
                                                directionality,
                                                complex_id=complex_id,
                                                spontaneous=spontaneous,
                                                update=update, keff=keff)
