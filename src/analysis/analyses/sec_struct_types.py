from src.analysis.analysis import Analysis
from src.rna_structure.structure_io import StructureIO
from src.rna_folding.rna_folder import RNAFolder


class SecondaryStructureTypes(Analysis):
    """
    Reads an input structure file (connectivity table or
    dot bracket) and returns the percent of bases in the
    following categories of secondary structure: watson-crick
    stem, non-WC stem, pseudoknot, or no stem.

    Parameters
    ----------
    config : AnalysisParser
        Object containing user inputs

    """

    def __init__(self, config):
        super().__init__(config)
        self.connect_table = StructureIO()._ct_to_dataframe(self.config.args.input)
        self.num_bases = self.connect_table["Index"].iloc[-1]
        self.wc_interactions = RNAFolder(config).interactions
        # initialize counters
        self.wc_base_pairs = 0
        self.nonwc_base_pairs = 0
        self.pseudoknots = 0
        self.no_pair = 0
        self._analyze()

    def _analyze(self):
        for i in range(self.num_bases):
            if self.connect_table["Paired With"].iloc[i] == 0:
                self.no_pair += 1
            elif self.connect_table["Paired With"].iloc[i] != 0:
                if self._check_wc_interactions(
                    self.connect_table["Nucleotide"].iloc[i],
                    self.connect_table["Nucleotide"].iloc[
                        self.connect_table["Paired With"].iloc[i] - 1
                    ],
                ):
                    self.wc_base_pairs += 1
                else:
                    self.nonwc_base_pairs += 1
        print(self.no_pair, self.wc_base_pairs, self.nonwc_base_pairs)

    def _check_wc_interactions(self, base1, base2):
        for pair in self.wc_interactions:
            if pair[0] == base1 and pair[1] == base2:
                return True
        return False
