from src.analysis.analysis import Analysis
from src.structure_io.structure_io import StructureIO


class SecondaryStructureTypes(Analysis):
    """
    Reads an input structure file (connectivity table or
    dot bracket) and returns the percent of bases in the
    following categories of secondary structure: watson-crick
    stem, non-WC stem, pseudoknot, or no stem.

    Parameters
    ----------
    config : AnalysisParser
        Object containing user inputs

    """

    def __init__(self, config):
        super().__init__(config)
        self.connect_table = StructureIO()._ct_to_dataframe(self.config.args.input)
        self.num_bases = self.connect_table["Index"].iloc[-1]
        # initialize counters
        self.wc_base_pairs = 0
        self.nonwc_base_pairs = 0
        self.pseudoknots = 0
        self.no_pair = 0

    def _analyze(self):
        pass
