from src.analysis.analysis import Analysis
from src.rna_folding.rna_folder import RNAFolder
from src.rna_structure.structure import RNAStructure
from src.rna_structure.structure_io import StructureIO
from src.rna_structure.structure_convert import StructureConvert


class ComputeEnergy(Analysis):
    """
    For a given input sequence and structure (connectivity table or TODO dot-bracket),
    calculate energy according to the Hamiltonian defined in Fox et al. PLoS One. 2022.
    https://doi.org/10.1371/journal.pcbi.1010032.

    Parameters
    ----------
    config : AnalysisParser
        Object containing user inputs

    """

    def __init__(self, config):
        super().__init__(config)
        # read sequence from structure file
        self.connect_table = StructureIO()._ct_to_dataframe(self.config.args.input)
        seq = StructureIO()._get_sequence_from_connect_table(self.connect_table)

        # generate stems for sequence - need min stem length and min loop length!
        self.rna_folder_obj = RNAFolder(config)
        self.rna_folder_obj._declare_stem_vars(seq)

        # convert from connectivity table to stem tuples
        self.rna_struct_obj = RNAStructure()
        self.struct_conv_obj = StructureConvert()
        self.observed_stems = self.struct_conv_obj._connect_table_to_stems(
            self.rna_folder_obj.n, self.connect_table
        )

        # k_neighbor_energy inherits from this analysis. below are actions
        # which are only relevant to computing the energy of a structure
        if self.config.args.command == "compute_energy":
            self.rna_folder_obj.stems = self.observed_stems
            self.rna_folder_obj.len_stem_list = len(self.observed_stems)
            self.active_stem_indices = [
                i for i in range(self.rna_folder_obj.len_stem_list)
            ]
            self.rna_folder_obj._compute_h_and_J()
            self._analyze()

    def _analyze(self):
        self.score = self.rna_folder_obj._calc_score(self.active_stem_indices)
        self.config.log.info("Sequence: " + self.rna_folder_obj.nseq)
        self.config.log.info("Sequence Length: " + str(self.rna_folder_obj.n))
        self.config.log.info("Energy of input structure: " + str(self.score))
        print(self.rna_folder_obj.n, self.score)
