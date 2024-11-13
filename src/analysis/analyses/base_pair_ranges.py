from src.analysis.analysis import Analysis
from src.rna_structure.structure_io import StructureIO
import numpy as np


class BasePairRanges(Analysis):
	"""
	Determines the average, minimum, and maximum base pair lengths,
	as well as sequence length for a connectivity file.

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
		interactions = np.abs(input_df.loc[paired_indices, "Paired With"].values - input_df.loc[paired_indices, "Index"].values)

		avg_range = np.mean(interactions)
		min_range = np.min(interactions)
		max_range = np.max(interactions)
		seq_len = input_df["Index"].iloc[-1]
		
		self.config.log.info('Average pair range, Min pair range, Max pair range, Sequence Length')

		outputs = (avg_range, min_range, max_range, seq_len)		
		vals = ", ".join([str(outputs[k]) for k in range(len(outputs))])
		self.config.log.info(vals)

		print(avg_range, min_range, max_range, seq_len)





