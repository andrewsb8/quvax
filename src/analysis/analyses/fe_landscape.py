from src.analysis.analysis import Analysis


class FreeEnergyLandscape(Analysis):
    """
    Characterizing the free energy local free energy landscape explored of an optimization with design.py

    Parameters
    ----------
    """
    def __init__(self, config):
        super().__init__(config)
        self._analyze()

    def _analyze(self):
        pass
