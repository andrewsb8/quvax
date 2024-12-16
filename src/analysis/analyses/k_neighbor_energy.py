from src.analysis.analyses.compute_energy import ComputeEnergy


class KNeighborEnergySearch(ComputeEnergy):
    """
    For a given input sequence and structure (connectivity table or TODO dot-bracket),
    calculate distribution of energies by changing k neighbors (stems)

    Parameters
    ----------
    config : AnalysisParser
        Object containing user inputs

    """

    def __init__(self, config):
        super().__init__(config)
        self.rna_folder_obj._gen_stems()

        # find indices of stems from connect table and place in sorted list
        self.active_stem_indices = self._find_observed_stem_indices(
            self.observed_stems, self.rna_folder_obj.stems
        )

        self.rna_folder_obj.len_stem_list = len(self.rna_folder_obj.stems)
        self.rna_folder_obj._compute_h_and_J()
        # keep list of calculated energies, starting with the initial structure
        if self.active_stem_indices != []:
            self.energies = [self.rna_folder_obj._calc_score(self.active_stem_indices)]
        else:
            self.energies = [self.rna_folder_obj.no_stem_penalty]
        self._analyze()

    def _analyze(self):
        # brute force solution for k = 1
        new_stems = []
        valid_neighbor_count = 0
        lower_energy_neighbor_count = 0
        for i in range(self.rna_folder_obj.len_stem_list):
            if i in self.active_stem_indices:
                new_stems = [j for j in self.active_stem_indices if j != i]
            elif i not in self.active_stem_indices:
                new_stems = self.active_stem_indices + [i]

            if new_stems != []:
                new_energy = self.rna_folder_obj._calc_score(new_stems)
            else:
                new_energy = self.rna_folder_obj.no_stem_penalty

            if new_energy < 0:
                valid_neighbor_count += 1
            if new_energy < self.energies[0]:
                lower_energy_neighbor_count += 1
            self.energies.append(new_energy)
        self.config.log.info("Sequence: " + self.rna_folder_obj.nseq)
        self.config.log.info("Sequence Length: " + str(self.rna_folder_obj.n))
        self.config.log.info(
            "Outputs: energy of input structure, count of neighbors with energy < 0, count of neighbors with lower energy than input structure, minimum energy of neighbors"
        )
        self.config.log.info(
            str(self.energies[0])
            + " "
            + str(valid_neighbor_count)
            + " "
            + str(lower_energy_neighbor_count)
            + " "
            + str(min(self.energies))
        )
        output_energies = ", ".join(str(_) for _ in self.energies)
        self.config.log.info("List of input structure and neighbor output energies")
        self.config.log.info(output_energies)
        print(output_energies)
