from src.analysis.analysis import Analysis
from src.rna_structure.structure_io import StructureIO
from src.rna_structure.structure_convert import StructureConvert
from src.rna_structure.structure import RNAStructure


class DBtoCT(Analysis):
    """
    Convert dot-bracket to connectivity table

    Parameters
    __________

    config : AnalysisParser
            Object containing user inputs

    """

    def __init__(self, config):
        super().__init__(config)
        self._analyze()

    def _analyze(self):
        sio = StructureIO()
        sc = StructureConvert()
        self.seq, self.dot_bracket = sio._read_dot_bracket(self.config.args.input)
