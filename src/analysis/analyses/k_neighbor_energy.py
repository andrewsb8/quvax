from src.analysis.analysis import Analysis
from src.rna_folding.rna_folder import RNAFolder
from src.rna_structure.structure import RNAStructure
from src.rna_structure.structure_convert import StructureConvert


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
        # read sequence and stems from structure file
        # generate stems for sequence
        # find stem indices and keep sorted list of them
        # calculate the hamiltonian matrices
        # keep list of calculated energies, starting with the initial structure
        self._analyze()

    def _analyze(self):
        # for now, write loop to change k (keep it <=2) stems (add or delete)
        # this should have checks to maybe prevent duplicate checks.... not sure what kind yet
        # calculate energy of changed structure and record it
        pass
