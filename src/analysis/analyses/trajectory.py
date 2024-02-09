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
    """

    def __init__(self, config):
        super().__init__(config)
        self._analyze()

    def _analyze(self):
        number_of_gens = int(
            len(self.config.data["sequences"]) / self.config.data["generation_size"]
        )
        self.iterations = np.arange(1, number_of_gens)

        for m in range(self.config.data["generation_size"]):
            self.energy_diff = [
                self.config.data["energies"][
                    m + (k * self.config.data["generation_size"])
                ]
                for k in range(number_of_gens)
            ]

            file = self.config.args.output + "-" + str(m)

            self._print_output_2D(file, [self.iterations, self.energy_diff])
