from __future__ import print_function, division, absolute_import

from warnings import warn
from six import iteritems

import pandas
from Bio import SeqIO

import cobrame
from cobrame.util import dogma


def add_transcription_reaction(me_model, tu_name, locus_ids, sequence,
                               update=True):
    """
    Create TranscriptionReaction object and add it to ME-Model.
    This includes the necessary transcription data.

    Parameters
    ----------
    me_model : :class:`cobrame.core.model.MEModel`
        The MEModel object to which the reaction will be added

    tu_name : str
        ID of TU being transcribed.
        The TranscriptionReaction will be added as "transcription_+TU_name"
        The TranscriptionData will be added as just "TU_name"

    locus_ids : set
        Set of locus IDs that the TU transcribes

    sequence : str
        Nucleotide sequence of the TU.

    update : bool
        If True, use TranscriptionReaction's update function to update and
        add reaction stoichiometry

    Returns
    -------
    :class:`cobrame.core.reaction.TranscriptionReaction`
        TranscriptionReaction for the TU
    """

    transcription = cobrame.TranscriptionReaction("transcription_" + tu_name)
    transcription.transcription_data = \
        cobrame.TranscriptionData(tu_name, me_model)
    transcription.transcription_data.nucleotide_sequence = sequence
    transcription.transcription_data.RNA_products = {"RNA_" + i
                                                     for i in locus_ids}

    me_model.add_reaction(transcription)
    if update:
        transcription.update()
    return transcription


def create_transcribed_gene(me_model, locus_id, rna_type, seq,
                            left_pos=None, right_pos=None, strand=None):
    """
     Creates a `TranscribedGene` metabolite object and adds it to the ME-model

    Parameters
    ----------
    me_model : :class:`cobrame.core.model.MEModel`
        The MEModel object to which the reaction will be added

    locus_id : str
        Locus ID of RNA product.
        The TranscribedGene will be added as "RNA + _ + locus_id"

    left_pos : int or None
        Left position of gene on the sequence of the (+) strain

    right_pos : int or None
        Right position of gene on the sequence of the (+) strain

    seq : str
        Nucleotide sequence of RNA product.
        Amino acid sequence, codon counts, etc. will be calculated based on
        this string.

    strand : str or None
        - (+) if the RNA product is on the leading strand
        - (-) if the RNA product is on the complementary strand

    rna_type : str
        Type of RNA of the product.
        tRNA, rRNA, or mRNA
        Used for determining how RNA product will be processed.

    Returns
    -------
        :class:`cobrame.core.component.TranscribedGene`
            Metabolite object for the RNA product
    """
    gene = cobrame.TranscribedGene('RNA_' + locus_id, rna_type, seq)
    gene.left_pos = left_pos
    gene.right_pos = right_pos
    gene.strand = strand

    me_model.add_metabolites([gene])
    return gene


def add_translation_reaction(me_model, locus_id, dna_sequence,
                             update=False):
    """
    Creates and adds a TranslationReaction to the ME-model as well as the
    associated TranslationData

    A dna_sequence is required in order to add a TranslationReaction to the
    ME-model

    Parameters
    ----------
    me_model : :class:`cobra.core.model.MEModel`
        The MEModel object to which the reaction will be added

    locus_id : str
        Locus ID of RNA product.
        The TranslationReaction will be added as "translation + _ + locus_id"
        The TranslationData will be added as "locus_id"

    dna_sequence : str
        DNA sequence of the RNA product. This string should be reverse
        transcribed if it originates on the complement strand.

    update : bool
        If True, use TranslationReaction's update function to update and
        add reaction stoichiometry

    """

    # Create TranslationData
    translation_data = \
        cobrame.TranslationData(locus_id, me_model, "RNA_" + locus_id,
                                "protein_" + locus_id)
    translation_data.nucleotide_sequence = dna_sequence

    # Add RNA to model if it doesn't exist
    if "RNA_" + locus_id not in me_model.metabolites:
        warn('RNA_%s not present in model. Adding it now.')
        rna = cobrame.TranscribedGene('RNA_' + locus_id, 'mRNA', dna_sequence)
        me_model.add_metabolites(rna)

    # Create and add TranslationReaction with TranslationData
    translation_reaction = \
        cobrame.TranslationReaction("translation_" + locus_id)
    me_model.add_reaction(translation_reaction)
    translation_reaction.translation_data = translation_data

    if update:
        translation_reaction.update()


def convert_aa_codes_and_add_charging(me_model, trna_aa, trna_to_codon,
                                      verbose=True):
    """
    Adds tRNA charging reactions for all tRNAs in ME-model

    Parameters
    ----------
    me_model : :class:`cobra.core.model.MEModel`
        The MEModel object to which the reaction will be added

    trna_aa : dict
        Dictionary of tRNA locus ID to 3 letter codes of the amino acid
        that the tRNA contributes

        {tRNA identifier (locus_id): amino_acid_3_letter_code}

    trna_to_codon : dict
        Dictionary of tRNA identifier to the codon which it associates

        {tRNA identifier (locus_id): codon_sequence}

    verbose : bool
        If True, display metabolites that were not previously added to the
        model and were thus added when creating charging reactions
    """
    # convert amino acid 3 letter codes to metabolites
    for tRNA, aa in list(iteritems(trna_aa)):
        if aa == "OTHER":
            trna_aa.pop(tRNA)
        elif aa == "Sec":
            # Charge with precursor to selenocysteine
            trna_aa[tRNA] = me_model.metabolites.get_by_id('cys__L_c')
        elif aa == "Gly":
            trna_aa[tRNA] = me_model.metabolites.get_by_id("gly_c")
        else:
            trna_aa[tRNA] = \
                me_model.metabolites.get_by_id(aa.lower() + "__L_c")

    # add in all the tRNA charging reactions
    for tRNA, aa in iteritems(trna_aa):
        for codon in trna_to_codon[tRNA]:
            trna_data = cobrame.tRNAData("tRNA_" + tRNA + "_" + codon,
                                         me_model, aa.id, "RNA_" + tRNA, codon)
            charging_reaction = \
                cobrame.tRNAChargingReaction("charging_tRNA_%s_%s" % (tRNA,
                                                                      codon))
            charging_reaction.tRNA_data = trna_data

            me_model.add_reaction(charging_reaction)
            charging_reaction.update(verbose=verbose)


def build_reactions_from_genbank(me_model, gb_filename, tu_frame=None,
                                 element_types={'CDS', 'rRNA',
                                                'tRNA', 'ncRNA'},
                                 verbose=True, frameshift_dict=None,
                                 trna_to_codon=None, update=True):

    # TODO handle special RNAse without type ('b3123')
    """Creates and adds transcription and translation reactions using genomic
     information from the organism's genbank file. Adds in the basic
     requirements for these reactions. Organism specific components are added
     ...

    Parameters
    ----------
    me_model : :class:`cobrame.core.model.MEModel`
        The MEModel object to which the reaction will be added

    gb_filename : str
        Local name of the genbank file that will be used for ME-model
        construction

    tu_frame : :class:`pandas.DataFrame`
        DataFrame with indexes of the transcription unit name and columns
        containing the transcription unit starting and stopping location on
        the genome and whether the transcription unit is found on the
        main (+) strand or complementary (-) strand.

        If no transcription unit DataFrame is passed into the function,
        transcription units are added corresponding to each transcribed
        gene in the genbank file.

    element_types : set
        Transcription reactions will be added to the ME-model for all RNA
        feature.types in this set. This uses the nomenclature of the
        genbank file (gb_filename)

    verbose : bool
        If True, display metabolites that were not previously added to the
        model and were thus added when creating charging reactions

    frameshift_dict : dict
        {locus_id: genome_position_of_TU}

        If a locus_id is in the frameshift_dict, update it's nucleotide
        sequence to account of the frameshift

    """
    if not frameshift_dict:
        frameshift_dict = {}

    if not trna_to_codon:
        trna_to_codon = {}

    metabolites = me_model.metabolites

    # Load genbank file and extract DNA sequence
    gb_file = SeqIO.read(gb_filename, 'gb')
    full_seq = str(gb_file.seq)

    # Dictionary of tRNA locus ID to the 3 letter code for the amino acid it
    # contributes
    trna_aa = {}

    # If no tu_frame is provided generate a new TU frame where each mRNA gets
    # its own TU
    using_tus = tu_frame is not None
    if not using_tus:
        tu_frame = pandas.DataFrame.from_dict(
            {"TU_" + i.qualifiers["locus_tag"][0]:
                {"start": int(i.location.start),
                 "stop": int(i.location.end),
                 "strand": "+" if i.strand == 1 else "-"}
             for i in gb_file.features if
             i.type in element_types},
            orient="index")

    # Create transcription reactions for each TU and DNA sequence.
    # RNA_products will be added so no need to update now
    for tu_id in tu_frame.index:
        # subtract 1 from TU start site to account for 0 indexing
        sequence = dogma.extract_sequence(full_seq, tu_frame.start[tu_id]-1,
                                          tu_frame.stop[tu_id],
                                          tu_frame.strand[tu_id])

        add_transcription_reaction(me_model, tu_id, set(), sequence,
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
        rna_type = 'mRNA' if feature.type == 'CDS' else feature.type
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
        gene = create_transcribed_gene(me_model, bnum, rna_type, seq, left_pos,
                                       right_pos, strand)

        # ---- Add translation reaction for mRNA ----
        if rna_type == "mRNA":
            add_translation_reaction(me_model, bnum, dna_sequence=seq)

        # ---- Create dict to use for adding tRNAChargingReactions ----
        # tRNA_aa = {'amino_acid':'tRNA'}
        elif rna_type == "tRNA":
            trna_aa[bnum] = feature.qualifiers["product"][0].split("-")[1]

        # ---- Associate TranscribedGene to a TU ----
        parent_tu = tu_frame[
            (tu_frame.start - 1 <= left_pos) & (tu_frame.stop >= right_pos) & (
             tu_frame.strand == strand)].index

        if len(parent_tu) == 0:
            if verbose:
                warn('No TU found for %s %s' % (rna_type, bnum))
            tu_id = "TU_" + bnum
            parent_tu = [tu_id]
            add_transcription_reaction(me_model, tu_id, set(), seq,
                                       update=False)

        for TU_id in parent_tu:
            me_model.process_data.get_by_id(TU_id).RNA_products.add(
                    "RNA_" + bnum)

    convert_aa_codes_and_add_charging(me_model, trna_aa, trna_to_codon,
                                      verbose=verbose)
    if update:
        # Update all newly added reactions
        for r in me_model.reactions:
            if isinstance(r, (cobrame.TranscriptionReaction,
                              cobrame.TranslationReaction)):
                r.update()


def add_m_model_content(me_model, m_model, complex_metabolite_ids=None):
    """
    Add metabolite and reaction attributes to me_model from m_model. Also
    creates StoichiometricData objects for each reaction in m_model, and adds
    reactions directly to me_model if they are exchanges or demands.

    Parameters
    ----------
    me_model : :class:`cobrame.core.model.MEModel`
        The MEModel object to which the content will be added

    m_model : :class:`cobra.core.model.Model`
        The m_model which will act as the source of metabolic content for
        MEModel

    complex_metabolite_ids : list
        List of complexes which are 'metabolites' in the m-model reaction
        matrix, but should be treated as complexes

    """
    if not complex_metabolite_ids:
        complex_metabolite_ids = []

    for met in m_model.metabolites:
        if met.id in complex_metabolite_ids:
            new_met = cobrame.Complex(met.id)
        elif met.id.startswith("RNA"):
            raise ValueError('Processed M-model should not contain RNAs (%s)' %
                             met.id)
        else:
            new_met = cobrame.Metabolite(met.id)
        new_met.name = met.name
        new_met.formula = met.formula
        new_met.compartment = met.compartment
        new_met.charge = met.charge
        new_met.annotation = met.annotation
        new_met.notes = met.notes
        me_model.add_metabolites(new_met)

    for reaction in m_model.reactions:
        if reaction.id.startswith("EX_") or reaction.id.startswith("DM_"):
            new_reaction = cobrame.MEReaction(reaction.id)
            me_model.add_reaction(new_reaction)
            new_reaction.lower_bound = reaction.lower_bound
            new_reaction.upper_bound = reaction.upper_bound
            for met, stoichiometry in iteritems(reaction.metabolites):
                new_reaction.add_metabolites(
                    {me_model.metabolites.get_by_id(met.id): stoichiometry})

        else:
            reaction_data = cobrame.StoichiometricData(reaction.id, me_model)
            reaction_data.lower_bound = reaction.lower_bound
            reaction_data.upper_bound = reaction.upper_bound
            reaction_data._stoichiometry = {k.id: v for k, v
                                            in iteritems(reaction.metabolites)}


def add_dummy_reactions(me_model, dna_seq, update=True):
    """
    Add all reactions necessary to produce a dummy reaction catalyzed by
    "CPLX_dummy".

    Parameters
    ----------
    me_model : :class:`cobrame.core.model.MEModel`
        The MEModel object to which the content will be added

    dna_seq : str
        DNA sequence of dummy gene. Should be representative of the average
        codon composition, amino acid composition, length of a gene in the
        organism being modeled

    update : bool
        If True, run update functions on all transcription, translation,
        complex formation, and metabolic reactions added when constructing
        dummy reactions.

    """
    dummy = cobrame.StoichiometricData("dummy_reaction", me_model)
    dummy.lower_bound = 0
    dummy.upper_bound = 1000
    dummy._stoichiometry = {'CPLX_dummy': -1}

    create_transcribed_gene(me_model, 'dummy', 'mRNA', dna_seq)
    add_transcription_reaction(me_model, "RNA_dummy", {"dummy"}, dna_seq)
    me_model.add_metabolites(cobrame.TranslatedGene("protein_" + "dummy"))
    add_translation_reaction(me_model, "dummy", dna_sequence=dna_seq,
                             update=update)
    try:
        complex_data = cobrame.ComplexData("CPLX_dummy", me_model)
    except ValueError:
        warn('CPLX_dummy already in model')
        complex_data = me_model.process_data.get_by_id('CPLX_dummy')
    complex_data.stoichiometry = {"protein_dummy": 1}
    if update:
        complex_data.create_complex_formation()


def add_complex_to_model(me_model, complex_id, complex_stoichiometry,
                         complex_modifications=None):
    """
    Adds ComplexData to the model for a given complex.

    Parameters
    ----------
    me_model : :class:`cobrame.core.model.MEModel`

    complex_id : str
        ID of the complex and thus the model ComplexData

    complex_stoichiometry : dict
        {complex_id: {protein_<locus_tag>: stoichiometry}}

    complex_modifications : dict
        {subreaction_id: stoichiometry}

    """

    if not complex_modifications:
        complex_modifications = {}

    complex_data = cobrame.ComplexData(complex_id, me_model)
    # must add update stoichiometry one by one since it is a defaultdict
    for metabolite, value in iteritems(complex_stoichiometry):
        complex_data.stoichiometry[metabolite] += value
    for modification, value in iteritems(complex_modifications):
        complex_data.subreactions[modification] = value


def add_subreaction_data(me_model, modification_id,
                         modification_stoichiometry,
                         modification_enzyme=None, verbose=True):
    """
    Creates a SubreactionData object for each modification defined by the
    function inputs.

    It's assumed every complex modification occurs spontaneously, unless a
    modification_enzyme argument is passed.

    If a modification uses an enzyme this can be updated after the
    SubreactionData object is already created


    Parameters
    ----------
        me_model : :class:`cobrame.core.model.MEModel`

    """

    if modification_id in me_model.process_data:
        if verbose:
            warn('Subreaction (%s) already in model' % modification_id)
        else:
            pass
    else:
        modification_data = cobrame.SubreactionData(modification_id, me_model)
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

    Parameters
    ----------
    me_model : :class:`cobrame.core.model.MEModel`

    complex_stoichiometry_dict : dict
        {unmodified_complex_id: {protein_<locus_tag>: stoichiometry}}

    complex_modification_dict : dict
        {modified_complex_id:{core_enzyme: unmodified_complex_id,
                             'modifications: {mod_metabolite: stoichiometry}}}

    """
    for complex_id, stoichiometry in iteritems(complex_stoichiometry_dict):
        add_complex_to_model(me_model, complex_id, stoichiometry, {})

    for modified_complex_id, info in iteritems(complex_modification_dict):
        modification_dict = {}
        for metabolite, number in iteritems(info['modifications']):
            modification_id = 'mod_' + metabolite

            # add modification as subreaction
            add_subreaction_data(me_model, modification_id, {metabolite: -1},
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

    Parameters
    ----------
    me_model : :class:`cobrame.core.model.MEModel`
        MEModel that the MetabolicReaction will be added to

    stoichiometric_data_id : str
        ID of the StoichiometricData for the reaction being added

    directionality : str
        - Forward: Add reaction that occurs in the forward direction
        - Reverse: Add reaction that occurs in the reverse direction

    complex_id : str or None
        ID of the ComplexData for the enzyme that catalyze the reaction
        being added.

    spontaneous : bool
        - If True and complex_id='' add reaction as spontaneous reaction
        - If False and complex_id='' add reaction as orphan (CPLX_dummy
          catalyzed)

    """
    # Get stoichiometric data for reaction being added
    try:
        stoichiometric_data = \
            me_model.process_data.get_by_id(stoichiometric_data_id)
    except KeyError:
        raise Exception("Stoichiometric data for %s has not been added to"
                        " model" % stoichiometric_data_id)

    # Get complex data and id based on arguments passed into function
    if type(complex_id) == str:
        complex_data = me_model.process_data.get_by_id(complex_id)
    elif complex_id is None and spontaneous is True:
        complex_id = "SPONT"
        complex_data = None
    elif complex_id is None and spontaneous is False:
        complex_id = "CPLX_dummy"
        try:
            complex_data = me_model.process_data.get_by_id(complex_id)
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

    r = cobrame.MetabolicReaction(''.join([stoichiometric_data_id,
                                          direction, complex_id]))
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

    Parameters
    ----------
    me_model : :class:`cobrame.core.model.MEModel`
        MEModel that the MetabolicReaction will be added to

    rxn_to_cplx_dict : dict
        {StoichiometricData.id: catalytic_enzyme_id}

    rxn_info_frame: :class:`pandas.Dataframe`
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
        complexes_list = rxn_to_cplx_dict.get(reaction_data.id, [None])

        # Add metabolic reactions for each isozyme
        for complex_id in complexes_list:
            directionality_list = []
            if reaction_data.lower_bound < 0:
                directionality_list.append('reverse')
            if reaction_data.upper_bound > 0:
                directionality_list.append('forward')
            elif reaction_data.upper_bound == 0 and \
                    reaction_data.lower_bound == 0:
                directionality_list.append('forward')
                warn('Reaction (%s) cannot carry flux' % reaction_data.id)
            for directionality in directionality_list:
                add_metabolic_reaction_to_model(me_model, reaction_data.id,
                                                directionality,
                                                complex_id=complex_id,
                                                spontaneous=spontaneous,
                                                update=update, keff=keff)
