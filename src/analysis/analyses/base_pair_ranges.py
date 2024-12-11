from src.analysis.analysis import Analysis
from src.rna_structure.structure_io import StructureIO
import numpy as np


class BasePairRanges(Analysis):
    """
    Determines the average, minimum, and maximum base pair lengths,
    in number of bases, for an input connectivity table.

    Parameters
    __________

    config : AnalysisParser
            Object containing user inputs

    """

    def __init__(self, config):
        super().__init__(config)
        self._analyze()

    def _analyze(self):
        input_df = StructureIO()._ct_to_dataframe(self.config.args.input)
        paired_indices = input_df["Paired With"] != 0
        interactions = np.abs(
            input_df.loc[paired_indices, "Paired With"].values
            - input_df.loc[paired_indices, "Index"].values
        )

        self.avg_range = np.mean(interactions)
        self.min_range = np.min(interactions)
        self.max_range = np.max(interactions)
        self.seq_len = len(input_df.index)

        self.config.log.info(
            "Outputs: Average pair range [# of bases], Min pair range, Max pair range, Sequence Length"
        )
        outputs = (self.avg_range, self.min_range, self.max_range, self.seq_len)
        vals = ", ".join([str(outputs[k]) for k in range(len(outputs))])
        self.config.log.info(vals)
        print(vals)
