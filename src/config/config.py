from abc import ABC, abstractmethod

class Config(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def _parse(self):
        pass
