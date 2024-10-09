from src.analysis.analysis import Analysis


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

    def _analyze(self):
        pass
