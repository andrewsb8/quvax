from src.analysis.analysis import Analysis
import difflib
import numpy as np


class FreeEnergyGeneration(Analysis):
    """
    Characterizing the free energy of mRNA sequences throughout the optimization
    process for a given protein sequence sampled with design.py.

    Parameters
    ----------
    """

    def __init__(self, config):
        super().__init__(config)
        self._analyze()

    def _analyze(self):
        self.codon_diff = [
            sum(
                [
                    self._calc_codon_diff(
                        self.config.data["sequences"][j][i * 3 : (i * 3) + 3],
                        self.config.data["sequences"][j+self.config.data["generation_size"]][i * 3 : (i * 3) + 3],
                    )
                    for i in range(int(len(self.config.data["sequences"][0]) / 3))
                ]
            )
            for j in range(len(self.config.data["sequences"])-self.config.data["generation_size"])
        ]
        self.energy_diff = [
            self._calc_energy_diff(energy) for energy in self.config.data["energies"]
        ]

        self._generate_output_2D([self.codon_diff, self.energy_diff])

    def _calc_codon_diff(self, mfe_codon, codon):
        """
        Checks if codons in two sequences are different.

        """
        diff = difflib.context_diff(mfe_codon, codon)
        for i, s in enumerate(diff):
            if s[0] == "+" or s[0] == "-":
                return 1
        return 0

    def _calc_energy_diff(self, energy):
        return energy - self.mfe
