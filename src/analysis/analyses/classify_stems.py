from src.analysis.analysis import Analysis
from src.rna_structure.structure_io import StructureIO
from src.rna_structure.structure_convert import StructureConvert
from src.rna_structure.structure import RNAStructure
import numpy as np
import warnings

class ClassifyStems(Analysis):
	"""
	Determines average stem length, minimum stem length, maximum stem length, 
	sequence length, total number of stems, number of stems in pseudoknots, and 
	number of overlapping stems for an input connectivity table.

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

		stems = StructureConvert()._connect_table_to_stems(sequence_len, connect_table_df)

		self.num_stems = len(stems)
		stem_lengths = list(inner_list[2] for inner_list in stems)
		self.min_stem = min(stem_lengths)
		self.max_stem = max(stem_lengths)
		self.seq_len = sequence_len
		self.avg_stem = sum(stem_lengths)/(self.num_stems)
		self.pseudos, self.overlaps = self._get_pseudos_and_overlap(stems)

		if self.overlaps > 0:
			warnings.warn("One or more stems overlap",category=UserWarning)
		self.config.log.info(
			"Outputs: Avg stem length, Min stem length, Max stem length, Sequence length, Number of stems, Number of pseudoknot stems, Number of overlapping Stems"
		)
		outputs = (self.avg_stem, self.min_stem, self.max_stem, self.seq_len, self.num_stems, self.pseudos, self.overlaps)
		vals = ", ".join([str(outputs[k]) for k in range(len(outputs))])
		self.config.log.info(vals)
		print(vals)


	def _get_pseudos_and_overlap(self,stems):
		pseudos = 0
		overlaps = 0
		for i in stems:
			for j in stems:
				if i != j:
					if RNAStructure()._is_pseudo(i,j)==True:
						pseudos += 1			
					if RNAStructure()._detect_stem_overlap(i,j)==True:
						overlaps += 1
		return pseudos, overlaps			
