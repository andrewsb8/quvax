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
        self.active_stem_indices = []
        for stem in self.observed_stems:
            found = False
            for j in range(len(self.rna_folder_obj.stems)):
                if stem == self.rna_folder_obj.stems[j]:
                    found = True
                    self.active_stem_indices.append(j)
            # if stem in connect table is not found in possible stem list
            # (self.rna_folder_obj.stems), likely because value for -ms is
            # higher than some stem lengths in input structure or a loop
            # smaller than specified by -ml is observed. Then add
            # observed stem to stem list and active stem list. This way,
            # the observed stems will be considered in the energy calculation
            # (or bit flipping process in k_neighbor_energy)
            if not found:
                self.rna_folder_obj.stems.append(stem)
                self.active_stem_indices.append(len(self.rna_folder_obj.stems) - 1)
        self.active_stem_indices.sort()
        print(self.rna_folder_obj.stems)
        print(self.active_stem_indices)

        self.rna_folder_obj.len_stem_list = len(self.rna_folder_obj.stems)
        self.rna_folder_obj._compute_h_and_J()
        # keep list of calculated energies, starting with the initial structure
        self.energies = [self.rna_folder_obj._calc_score(self.active_stem_indices)]
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

            new_energy = self.rna_folder_obj._calc_score(new_stems)
            if new_energy < 0:
                valid_neighbor_count += 1
            if new_energy < self.energies[0]:
                lower_energy_neighbor_count += 1
            self.energies.append(new_energy)
        # print(self.energies)
        print(
            self.energies[0],
            valid_neighbor_count,
            lower_energy_neighbor_count,
            min(self.energies),
        )
