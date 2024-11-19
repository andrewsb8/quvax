from src.analysis.analysis import Analysis
from src.rna_structure.structure_io import StructureIO
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
		stems = []
		i=0

		while i < sequence_len:
	            # move on from unpaired bases and don't double count base pairs
			if (
				connect_table_df["Paired With"].iloc[i]
				< connect_table_df["Index"].iloc[i]
				or connect_table_df["Paired With"].iloc[i] == 0
			):
				i += 1
                	# if base pair is found, need to define the length of stem  
			elif connect_table_df["Paired With"].iloc[i] != 0:
                	# loop through connect table until nonsequential base pair is
                	# found and then append to stems
				for j in range(i + 1, sequence_len - 1):
                    	# if nonsequential base pair ordering or a unpaired base is found, use information to generate stem
					if (
						connect_table_df["Paired With"].iloc[j]
						!= connect_table_df["Paired With"].iloc[j - 1] - 1
						or connect_table_df["Paired With"].iloc[j] == 0
					):
						stems.append(
							(
								connect_table_df["Index"].iloc[i],
								connect_table_df["Paired With"].iloc[i],
								j - i,
							)
						)
						i += j-i
						break

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



