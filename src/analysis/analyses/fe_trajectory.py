from src.analysis.analysis import Analysis
import numpy as np


class FETrajectory(Analysis):
    """
    Plotting the folding free energies as a function of codon optimization step.
    This analysis will produce a number of output files equal to the population
    size (n_trials in design.py).

    Parameters
    ----------
    config : AnalysisParser
        Object containing user inputs

    """

    def __init__(self, config):
        super().__init__(config)
        self._analyze()

    def _analyze(self):
        self.iterations = np.arange(
            0, self.config.sim_details["generations_sampled"] + 1
        )

        for m in range(self.config.sim_details["generation_size"]):
            self.config.db_cursor.execute(
                f"SELECT energies from OUTPUTS WHERE population_key = {m};"
            )
            self.data = self.config.db_cursor.fetchall()
            self.energies = [dat[0] for dat in self.data]

            file = self.config.args.output + "-" + str(m)
            self._print_output_2D(file, [self.iterations, self.energies])
