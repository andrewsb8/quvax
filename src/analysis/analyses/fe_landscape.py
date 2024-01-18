from src.analysis.analysis import Analysis


class FreeEnergyLandscape(Analysis):
    """
    Characterizing the local free energy landscape of mRNA sequence space for a given protein sequence sampled with design.py.

    Parameters
    ----------
    """
    def __init__(self, config):
        super().__init__(config)
        self._analyze()

    def _analyze(self):
        pass
