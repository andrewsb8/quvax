from src.analysis.analysis import Analysis
from src.rna_folding.rna_folder import RNAFolder
from src.rna_structure.structure import RNAStructure
from src.rna_structure.structure_io import StructureIO
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
        self.connect_table = StructureIO()._ct_to_dataframe(self.config.args.input)
        seq = "".join([self.connect_table["Nucleotide"].iloc[i] for i in range(self.connect_table["Index"].iloc[-1])])

        # generate stems for sequence - need min stem lenght and min loop length!
        self.rna_folder_obj = RNAFolder(config) # no config needed here
        self.rna_folder_obj._fold_prep(seq)

        # convert from connectivity table to dot bracket for easy identification of pseudoknots
        self.rna_struct_obj = RNAStructure()
        self.struct_conv_obj = StructureConvert()
        stems = self.struct_conv_obj._connect_table_to_stems(
            self.rna_folder_obj.n, self.connect_table
        )

        # find stem indices and keep sorted list of them
        # keep list of calculated energies, starting with the initial structure
        self._analyze()

    def _analyze(self):
        # for now, write loop to change k (keep it <=2) stems (add or delete)
        # this should have checks to maybe prevent duplicate checks.... not sure what kind yet
        # calculate energy of changed structure and record it
        pass
