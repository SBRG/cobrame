from cobra import Metabolite as Component


class Metabolite(Component):
    pass


class TranscribedGene(Component):
    pass


class TranslatedGene(Component):
    pass


class Complex(Component):
    pass


class Ribosome(Complex):
    pass


class GenerictRNA(Component):
    pass


class RNAP(Complex):
    pass


class Constraint(Component):
    pass


class ModComplex(Component):
    pass


def create_component(component_id, default_type=Component):
    """creates a component and attempts to set the correct type"""
    if not isinstance(component_id, str):
        raise TypeError("%s must be a str, not %s" %
                        (repr(component_id), str(type(component_id))))
    if component_id.startswith("protein"):
        return TranslatedGene(component_id)
    elif component_id.startswith("RNA"):
        return TranscribedGene(component_id)
    elif component_id.startswith("ribosome"):
        return Ribosome(component_id)
    elif component_id.startswith("RNA_Polymerase"):
        return RNAP(component_id)
    else:
        return default_type(component_id)
