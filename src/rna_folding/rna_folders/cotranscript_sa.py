from src.rna_folding.rna_folders.simulated_annealer import SimulatedAnnealer
import math


class CoTranscriptSA(SimulatedAnnealer):
    def __init__(self, config):
        super().__init__(config)

    def _fold(self, sequence, post_process=False):
        transcript_rate = 10
        old_stems = []
        self.stems_used = []

        for i in range(math.ceil(len(sequence) / transcript_rate)):
            if (i + 1) * transcript_rate >= len(sequence):
                self._fold_prep(sequence)
            else:
                self._fold_prep(sequence[0 : (i + 1) * transcript_rate])
            if self.len_stem_list > 0:
                if i == 0:
                    self._compute_dwave_sa(nr=3, isg="random")
                if i > 0:
                    # add zeros (stem not used) to intial state for all
                    # new stems found from self._gen_stems() after extending
                    # the sequence
                    for j in range(len(old_stems), len(self.stems)):
                        self.stem_dict[j] = 0
                    self._compute_dwave_sa(nr=3, istates=self.stem_dict, isg="tile")
                self._stems_to_dot_bracket(self.n, self.stems_used)
            else:
                self.stem_dict = {}
                self._stems_to_dot_bracket(self.n, [])
            old_stems = self.stems

        if post_process:  # equates to "if call was from fold.py"
            self._stems_to_connect_list(self.n, self.stems_used)
            self._post_process()
