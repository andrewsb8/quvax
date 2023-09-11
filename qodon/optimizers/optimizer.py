from abc import ABC, abstractmethod
from include.parser import Parser

class Optimizer(ABC):
    """
    Parent class for all codon optimizer classes.

    Parameters
    ----------
    config : Parser
        Object containing user inputs

    """
    def __init__(self, config: Parser):
        self.config = config

    @abstractmethod
    def _optimize(self):
        pass
