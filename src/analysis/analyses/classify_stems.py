from src.analysis.analysis import Analysis
from src.rna_structure.structure_io import StructureIO
from src.rna_structure.structure_convert import StructureConvert
from src.rna_structure.structure import RNAStructure
import numpy as np


class ClassifyStems(Analysis):
    """
    Determines average stem length, minimum stem length, maximum stem length,
    sequence length, total number of stems, number of pseudoknots stems,
    and number of overlapping stems for an input connectivity table.

    Parameters
    __________

    config : AnalysisParser
            Object containing user inputs

    """

    def __init__(self, config):
        super().__init__(config)
        self._analyze()

    def _analyze(self):
        connect_table_df = StructureIO()._ct_to_dataframe(self.config.args.input)
        sequence_len = len(connect_table_df["Index"])

        self.stems = StructureConvert()._connect_table_to_stems(
            sequence_len, connect_table_df
        )

        self.rna_struct_obj = RNAStructure()

        self.num_stems = len(self.stems)
        stem_lengths = list(inner_list[2] for inner_list in self.stems)
        self.min_stem = min(stem_lengths)
        self.max_stem = max(stem_lengths)
        self.seq_len = sequence_len
        self.avg_stem = sum(stem_lengths) / (self.num_stems)
        self.pseudos, self.overlaps, self.stems_of_loops = self._get_pseudos_and_overlaps()

        self.config.log.info(
            "Outputs: Avg stem length, Min stem length, Max stem length, Sequence length, Number of stems, Number of pseudoknot stems, Number of overlapping stems, Number of stems formed by two hairpin loops"
        )
        outputs = (
            self.avg_stem,
            self.min_stem,
            self.max_stem,
            self.seq_len,
            self.num_stems,
            self.pseudos,
            self.overlaps,
            self.stems_of_loops,
        )
        vals = ", ".join([str(outputs[k]) for k in range(len(outputs))])
        self.config.log.info(vals)
        if self.overlaps > 0:
            self.config.log.warning("Overlap detected in one or more stems")
        print(vals)

    def _get_pseudos_and_overlaps(self):
        pseudos = 0
        overlaps = 0
        stems_of_loops = 0
        for i in range(self.num_stems):
            if self.rna_struct_obj._detect_stem_of_loops(self.stems[i], self.stems):
                stems_of_loops += 1
            for j in range(i + 1, self.num_stems):
                self.rna_struct_obj = RNAStructure()
                if self.rna_struct_obj._is_pseudo(self.stems[i], self.stems[j]):
                    pseudos += 1
                if self.rna_struct_obj._detect_stem_overlap(
                    self.stems[i], self.stems[j]
                ):
                    overlaps += 1
        return pseudos, overlaps, stems_of_loops
