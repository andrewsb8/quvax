from abc import ABC, abstractmethod


class Analysis(ABC):
    """
    Parent class for all analyses. Will read output file from optimization process.

    Parameters
    ----------
    """
    def __init__(self):
        pass

    @abstractmethod
    def _analyze(self):
        pass
