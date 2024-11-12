from src.analysis.analysis import Analysis
from src.rna_structure.structure_io import StructureIO
import numpy as np


class PairingRange(Analysis):
	"""
	Determines the average base pair length for a connectivity file.
	Outputs average length as a fraction of the sequence length. 

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
		scaled_avg = avg_range/(input_df["Index"].iloc[-1]+1)
		self.config.log.info(scaled_avg)
		print(scaled_avg)





