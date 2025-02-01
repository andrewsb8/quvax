from src.analysis.analyses.compute_energy import ComputeEnergy


class Unfold(ComputeEnergy):
    """
    Compute the energy of a structure then repeat by removing one stem at a time
    from the structure.

    Parameters
    ----------
    config : AnalysisParser
        Object containing user inputs

    """

    def __init__(self, config):
        super().__init__(config)

        if self.config.args.target_stem_length == -1:
            self.rna_folder_obj._gen_stems()
            self.active_stem_indices = self._find_observed_stem_indices(
                self.observed_stems, self.rna_folder_obj.stems
            )
            self.rna_folder_obj.len_stem_list = len(self.rna_folder_obj.stems)
        else:
            self.rna_folder_obj.stems = self.observed_stems
            self.rna_folder_obj.len_stem_list = len(self.observed_stems)
            self.active_stem_indices = [
                i for i in range(self.rna_folder_obj.len_stem_list)
            ]
        self.rna_folder_obj._compute_h_and_J()
        self._analyze()

    def _analyze(self):
        # brute force solution for k = 1
        stems = self.active_stem_indices
        energy = self.rna_folder_obj._calc_score(self.active_stem_indices)
        for i in range(len(self.observed_stems) - 1):
            print(len(stems), energy)
            stems.pop()
            energy = self.rna_folder_obj._calc_score(stems)
        print(len(stems), energy)
