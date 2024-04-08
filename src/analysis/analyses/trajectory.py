from src.analysis.analysis import Analysis
import numpy as np


class Trajectory(Analysis):
    """
    Plotting the free energy changes as a function of codon optimization step
    relative to the initial sequence of the population. This analysis will
    produce a number of output files equal to the population size (n_trials in
    design.py).

    Parameters
    ----------
    config : AnalysisParser
        Object containing user inputs

    """

    def __init__(self, config):
        super().__init__(config)
        self._analyze()

    def _analyze(self):
        self.iterations = np.arange(0, self.config.sim_details["number_generations"])

        for m in range(self.config.sim_details["number_generations"]):
            self.config.db_cursor.execute(f"SELECT sequences, energies from OUTPUTS WHERE population_key = {m};")
            data = self.config.db_cursor.fetchall()

            self.energy_diff = []

            """
            self.energy_diff = [
                self.config.data["energies"][
                    m + (k * self.config.data["generation_size"])
                ]
                for k in range(len(data))
            ]

            file = self.config.args.output + "-" + str(m)
            self._print_output_2D(file, [self.iterations, self.energy_diff])
            """
