from abc import ABC, abstractmethod

class Optimizer(ABC):
    """
    Parent class for all codon optimizer classes

    """
    @abstractmethod
    def __init__(self):
        pass
