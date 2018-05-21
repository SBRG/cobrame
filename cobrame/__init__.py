from __future__ import absolute_import

from distutils.version import StrictVersion
from warnings import warn

import cobra

from cobrame.core.component import (
    Metabolite, TranscribedGene, Constraint, TranslatedGene, ProcessedProtein,
    Complex, Ribosome, GenericComponent, GenerictRNA, RNAP)
from cobrame.core.model import MEModel
from cobrame.core.reaction import (
    MetabolicReaction, tRNAChargingReaction, ComplexFormation,
    PostTranslationReaction, TranscriptionReaction, TranslationReaction,
    GenericFormationReaction, SummaryVariable, MEReaction)
from cobrame.core.processdata import (
    StoichiometricData, SubreactionData, ComplexData, TranscriptionData,
    tRNAData, TranslationData, TranslocationData, GenericData,
    PostTranslationData)
from cobrame.util import mu


__version__ = "0.0.9"

if StrictVersion(cobra.__version__) > StrictVersion('0.5.11'):
    raise UserWarning('COBRApy version (%s) is not <= 0.5.11. We recommend '
                      'using verision 0.5.11' % cobra.__version__)
elif StrictVersion(cobra.__version__) < StrictVersion('0.5.11'):
    warn('COBRApy version is %s. We recommend using 0.5.11. Using earlier '
         'versions may cause errors' % cobra.__version__)
