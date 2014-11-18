from cobra import Model, DictList


class MEmodel(Model):
    def __init__(self, *args):
        Model.__init__(self, *args)
        self.metabolic_reaction_data = DictList()
        self.complex_data = DictList()
        self.translation_data = DictList()
        self.transcription_data = DictList()
