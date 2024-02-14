from src.analysis.analysis import Analysis
import numpy as np


class FreeEnergyLandscape(Analysis):
    """
    Characterizing the local free energy landscape, relative to the lowest
    energy sequence sampled, of mRNA sequence space for a given protein
    sequence sampled with design.py.

    Parameters
    ----------
    """

    def __init__(self, config):
        super().__init__(config)
        self._analyze()

    def _analyze(self):
        self.mfe = np.min(self.config.data["energies"])
        self.mfe_index = np.argmin(self.config.data["energies"])
        self.mfe_codons = self.config.data["sequences"][self.mfe_index]

        self.codon_diff = [
            sum(
                [
                    self._calc_codon_diff(
                        self.mfe_codons[i * 3 : (i * 3) + 3],
                        self.config.data["sequences"][j][i * 3 : (i * 3) + 3],
                    )
                    for i in range(int(len(self.mfe_codons) / 3))
                ]
            )
            for j in range(len(self.config.data["sequences"]))
        ]
        self.energy_diff = [
            self._calc_energy_diff(self.mfe, energy)
            for energy in self.config.data["energies"]
        ]

        self._print_output_2D(
            self.config.args.output, [self.codon_diff, self.energy_diff]
        )
