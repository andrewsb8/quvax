from src.analysis.analysis import Analysis
import numpy as np


class FreeEnergyLandscape(Analysis):
    """
    Characterizing the local free energy landscape, relative to the lowest
    energy sequence sampled, of mRNA sequence space for a given protein
    sequence sampled with design.py.

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

        self.config.db_cursor.execute(f"SELECT sequences, energies from OUTPUTS where sim_key = {self.config.sim_details['sim_key']};")
        data = self.config.db_cursor.fetchall()
        self.sequences = [dat[0] for dat in data]
        self.energies = [dat[1] for dat in data]
        self.mfe_sequence = self.sequences[np.argmin(self.energies)]

        self.codon_diff = self._codon_diff_list(self.sequences, self.mfe_sequence)
        self.energy_diff = [
            self._calc_energy_diff(self.config.sim_details["min_free_energy"], energy)
            for energy in self.energies
        ]

        self._print_output_2D(
            self.config.args.output, [self.codon_diff, self.energy_diff]
        )
