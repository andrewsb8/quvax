from src.analysis.analysis import Analysis
from src.rna_structure.structure import RNAStructure
from src.rna_structure.structure_io import StructureIO
from src.rna_structure.structure_convert import StructureConvert


class BasePairTypes(Analysis):
    """
    Reads an input structure file (connectivity table [TODO: and dot-bracket])
    and returns the number of bases in the following categories of secondary
    structure: watson-crick pairs, wobble pairs, non-WC pairs, pairs in
    pseudoknot, and unpaired bases.

    Parameters
    ----------
    config : AnalysisParser
        Object containing user inputs

    """

    def __init__(self, config):
        super().__init__(config)
        self.connect_table = StructureIO()._ct_to_dataframe(self.config.args.input)
        self.num_bases = len(self.connect_table.index)

        # convert from connectivity table to dot bracket for easy identification of pseudoknots
        self.rna_struct_obj = RNAStructure()
        self.struct_conv_obj = StructureConvert()
        stems = self.struct_conv_obj._connect_table_to_stems(
            self.num_bases, self.connect_table
        )
        self.dot_bracket = self.struct_conv_obj._stems_to_dot_bracket(
            self.num_bases, stems
        )
        self.config.log.info("Converted connectivity table to dot-bracket:")
        self.config.log.info(self.dot_bracket)
        self.wc_interactions = self.rna_struct_obj._get_wc_interactions()
        self.wobble_interactions = self.rna_struct_obj._get_wobble_interactions()

        # initialize counters
        self.wc_base_pairs = 0
        self.wobble_base_pairs = 0
        self.nonwc_base_pairs = 0
        self.pseudoknots = 0
        self.no_pair = 0
        self._analyze()

    def _analyze(self):
        for i in range(self.num_bases):
            if self.connect_table["Paired With"].iloc[i] == 0:
                self.no_pair += 1
            elif self.connect_table["Paired With"].iloc[i] != 0:
                if self._check_interactions(
                    self.wc_interactions,
                    self.connect_table["Nucleotide"].iloc[i],
                    self.connect_table["Nucleotide"].iloc[
                        self.connect_table["Paired With"].iloc[i] - 1
                    ],
                ):
                    self.wc_base_pairs += 1
                elif self._check_interactions(
                    self.wobble_interactions,
                    self.connect_table["Nucleotide"].iloc[i],
                    self.connect_table["Nucleotide"].iloc[
                        self.connect_table["Paired With"].iloc[i] - 1
                    ],
                ):
                    self.wobble_base_pairs += 1
                else:
                    self.nonwc_base_pairs += 1
                # check for pseudoknots only if there is a pair
                if self.dot_bracket[i] == "[" or self.dot_bracket[i] == "]":
                    self.pseudoknots += 1

        # print output to log and stdout
        self.config.log.info(
            "# unpaired bases, bases in Watson-Crick pairs, bases in Wobble pairs, bases in non-WC pairs, bases in pseudoknots, number of bases (sequence length)"
        )
        output = ", ".join(
            str(item)
            for item in [
                self.no_pair,
                self.wc_base_pairs,
                self.wobble_base_pairs,
                self.nonwc_base_pairs,
                self.pseudoknots,
                self.num_bases,
            ]
        )
        self.config.log.debug(output)
        print(output)

    def _check_interactions(self, interaction_list, base1, base2):
        for pair in interaction_list:
            if pair[0] == base1 and pair[1] == base2:
                return True
        return False
