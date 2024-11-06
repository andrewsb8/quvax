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
        # read sequence from structure file
        self.connect_table = StructureIO()._ct_to_dataframe(self.config.args.input)
        seq = "".join([self.connect_table["Nucleotide"].iloc[i] for i in range(self.connect_table["Index"].iloc[-1])])

        # generate stems for sequence - need min stem length and min loop length!
        self.rna_folder_obj = RNAFolder(config)
        self.rna_folder_obj._fold_prep(seq)

        # convert from connectivity table to stem tuples
        self.rna_struct_obj = RNAStructure()
        self.struct_conv_obj = StructureConvert()
        stems = self.struct_conv_obj._connect_table_to_stems(
            self.rna_folder_obj.n, self.connect_table
        )

        # find indices of stems from connect table and place in sorted list
        self.active_stem_indices = []
        for stem in stems:
            for j in range(len(self.rna_folder_obj.stems)):
                if stem == self.rna_folder_obj.stems[j]:
                    self.active_stem_indices.append(j)
        self.active_stem_indices.sort()

        # keep list of calculated energies, starting with the initial structure
        self.energies = [self.rna_folder_obj._calc_score(self.active_stem_indices)]

        self._analyze()

    def _analyze(self):
        # how to loop over every combination by brute force?
        # start with an initial condition, need to enumerate each possible combination
        pass
