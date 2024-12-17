from src.analysis.analysis import Analysis
from src.rna_structure.structure_io import StructureIO


class ContactOrder(Analysis):
    """
    Calculates the contact order of an input RNA structure
    from a connectivity table (TO DO: add dot-bracket) as
    follows:

        CO = ( 1/( Sequence Length * Number of Base Pairs) ) *
            Sum over base pairs i[
                Sum over base pairs j>i[ Distance between base pairs ]
            ]

    If no base pairs are in the input structure, then a value of inf is
    assigned.

    Parameters
    __________

    config : AnalysisParser
            Object containing user inputs

    """

    def __init__(self, config):
        super().__init__(config)

        self.connect_table = StructureIO()._ct_to_dataframe(self.config.args.input)
        self.seq = StructureIO()._get_sequence_from_connect_table(self.connect_table)
        self.seq_length = len(self.seq)
        self._analyze()

    def _analyze(self):
        base_pairs = self.connect_table.loc[
            (self.connect_table["Paired With"] != 0)
            & (self.connect_table["Paired With"] > self.connect_table["Index"])
        ]
        count_base_pairs = len(base_pairs.index)
        if count_base_pairs == 0:
            self.config.log.warning("No base pairs in input structure!")
            print(self.seq_length, "inf")
        else:
            self.contact_order = sum(
                base_pairs.loc[:, "Paired With"].to_numpy()
                - base_pairs.loc[:, "Index"].to_numpy()
            ) / (self.seq_length * count_base_pairs)
            print(self.seq_length, self.contact_order)
