from abc import ABC, abstractmethod

class Optimizer(ABC):
    """
    Parent class for all codon optimizer classes

    """
    def __init__(self, config):
        self.config = config

    @abstractmethod
    def _optimize(self):
        pass
