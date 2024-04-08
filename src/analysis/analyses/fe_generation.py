from src.analysis.analysis import Analysis
import numpy as np


class FreeEnergyGeneration(Analysis):
    """
    Characterizing the free energy of mRNA sequences throughout the optimization
    process for a given protein sequence sampled with design.py.

    Parameters
    ----------
    config : AnalysisParser
        Object containing user inputs

    """

    def __init__(self, config):
        super().__init__(config)
        self._analyze()

    def _analyze(self):
        self.codon_diff_all = []
        self.energy_diff_all = []

        for m in range(self.config.sim_details["generation_size"]):
            self.config.db_cursor.execute(f"SELECT sequences, energies from OUTPUTS WHERE population_key = {m};")
            data = self.config.db_cursor.fetchall()
            self.sequences = [dat[0] for dat in data]
            self.energies = [dat[1] for dat in data]

            self.codon_diff = [ sum(
                [
                    self._calc_codon_diff(
                        self.sequences[j][i * 3 : (i * 3) + 3],
                        self.sequences[j+1][i * 3 : (i * 3) + 3],
                    )
                    for i in range(int(len(self.sequences[0]) / 3))
                ]
            ) for j in range(len(self.sequences) - 1)]
            self.energy_diff = [self._calc_energy_diff(self.energies[i], self.energies[i+1]) for i in range(len(self.energies) - 1)]

            self.codon_diff_all.extend(self.codon_diff)
            self.energy_diff_all.extend(self.energy_diff)

        self._print_output_2D(
            self.config.args.output, [self.codon_diff_all, self.energy_diff_all]
        )
