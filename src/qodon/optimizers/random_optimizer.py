from src.qodon.optimizer import CodonOptimizer
import random
from operator import itemgetter
import numpy as np
import math


class RandomOptimizer(CodonOptimizer):
    """
    A basic implementation of a random optimizer for a codon sequence.

    Parameters
    ----------

    """
    def __init__(self, config):
        super().__init__(config)
        self._optimize()

    def _optimize(self):
        '''
        Main method for codon optimization

        '''

        num_extra_sequences = (self.config.args.codon_iterations * self.config.args.n_trials) - self.config.args.n_trials
        extra_sequences = self._generate_sequences(num_extra_sequences)
        self.initial_sequences.extend(extra_sequences)

        scores = [self._tf_fold(s) for s in self.initial_sequences]

        self.final_energies = scores
        self.mfe = np.min(scores)
        self.mfe_index = np.argmin(scores)

        self.final_codons = self._reverse_translate(self.initial_sequences)

        self._verify_dna(self.final_codons[self.mfe_index])

        print(self.mfe)
        print(self.final_codons[self.mfe_index])
