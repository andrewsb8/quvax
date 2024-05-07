from abc import ABC, abstractmethod
from src.params.design_parser import DesignParser


class RNAFolder(ABC):
    """
    Parent class for all RNA folding classes.

    Parameters
    ----------


    """

    def __init__(self, config):
        self.config = config
        self.interactions = [
            ("A", "U"),
            ("U", "A"),
            ("G", "C"),
            ("C", "G"),
            ("G", "U"),
            ("U", "G"),
        ]
        self.twobody_penalty = 500000
        self.pseudo_factor = 0.5
        self.no_stem_penalty = 500000
