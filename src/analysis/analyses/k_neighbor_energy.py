from src.analysis.analysis import Analysis


class kNeighborEnergySearch(Analysis):
    """
    For a given input sequence and structure (connectivity table or TODO dot-bracket),
    calculate distribution of energies by changing k neighbors (stems)

    Parameters
    ----------
    config : AnalysisParser
        Object containing user inputs

    """

    def __init__(self, config):
        super().__init__(config)
        self._analyze()

    def _analyze(self):
        pass
