from src.analysis.analysis import Analysis
from src.rna_structure.structure_io import StructureIO
from src.rna_structure.structure_convert import StructureConvert
import numpy as np

class StemLengths(Analysis):
	"""
	Determines the average stem length, minimum stem length, maximum stem length, 
	sequence length and total number of stems for an input connectivity table.

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
		connect_table_df.columns =["Index","Nucleotide","Previous","Next","Paired With","Counter"]
		sequence_len = len(connect_table_df["Index"])        	

		stems = StructureConvert()._connect_table_to_stems(sequence_len, connect_table_df)

		self.num_stems = len(stems)
		self.min_stem = min(inner_list[2] for inner_list in stems)
		self.max_stem = max(inner_list[2] for inner_list in stems)
		self.seq_len = sequence_len
		tot_length = sum(inner_list[2] for inner_list in stems)
		self.avg_stem = tot_length/(self.num_stems)
                                
		self.config.log.info(
			"Outputs: Avg stem length, Min stem length, Max stem length, Sequence Length, Number of Stems"
		)
		outputs = (self.avg_stem, self.min_stem, self.max_stem, self.seq_len, self.num_stems)
		vals = ", ".join([str(outputs[k]) for k in range(len(outputs))])
		self.config.log.info(vals)
		print(vals)



