from src.analysis.analysis import Analysis
import numpy as np


class CodonTrajectory(Analysis):
    """
    Plotting the number of differences in codon sequences as a function of codon
    optimization step. This analysis will produce a number of output files equal
    to the population size (n_trials in design.py).

    Parameters
    ----------
    config : AnalysisParser
        Object containing user inputs

    """

    def __init__(self, config):
        super().__init__(config)
        self._analyze()

    def _analyze(self):
        self.iterations = np.arange(0, self.config.sim_details["generations_sampled"])

        for m in range(self.config.sim_details["generation_size"]):
            self.config.db_cursor.execute(
                f"SELECT sequences from OUTPUTS WHERE population_key = {m} and sim_key = {self.config.sim_details['sim_key']};"
            )
            self.data = self.config.db_cursor.fetchall()
            self.sequences = [dat[0] for dat in self.data]

            self.codon_diff = self._codon_diff_list(self.sequences)

            file = self.config.args.output + "-" + str(m)
            self._print_output_2D(file, [self.iterations, self.codon_diff])
