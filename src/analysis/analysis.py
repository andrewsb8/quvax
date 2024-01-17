from abc import ABC, abstractmethod
from src.params.analysis_parser import AnalysisParser


class Analysis(ABC):
    """
    Parent class for all analyses. Will read output file from optimization process.

    Parameters
    ----------
    """
    def __init__(self, config: AnalysisParser):
        pass

    @abstractmethod
    def _analyze(self):
        pass
