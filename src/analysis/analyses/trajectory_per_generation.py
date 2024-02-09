from src.analysis.analysis import Analysis
import numpy as np


class TrajectoryPerGeneration(Analysis):
    """
    Plotting the free energy changes as a function of codon optimization step
    in subsequent generations. This analysis will produce a number of output
    files equal to the population size (n_trials in design.py).

    Parameters
    ----------
    """

    def __init__(self, config):
        super().__init__(config)
        self._analyze()

    def _analyze(self):
        self.iterations = np.arange(
            1, len(self.config.data["sequences"]) / self.config.data["generation_size"]
        )

        self.energy_diff = [
            self._calc_energy_diff(
                self.config.data["energies"][k],
                self.config.data["energies"][k + self.config.data["generation_size"]],
            )
            for k in range(
                0,
                len(self.config.data["sequences"])
                - self.config.data["generation_size"],
                self.config.data["generation_size"],
            )
        ]

        self._print_output_2D(
            self.config.args.output, [self.iterations, self.energy_diff]
        )

    def _calc_energy_diff(self, ref_energy, energy):
        return energy - ref_energy
