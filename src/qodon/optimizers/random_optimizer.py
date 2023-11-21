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
        tmp = []
        for i in range(len(self.initial_sequences)):
            tmp.append([self._reverse_translate(self.initial_sequences[i])])

        scores = [self._fold_rna(s) for s in tmp]

        self._extend_output(self.initial_sequences, scores, None)
        self._get_optimized_sequence()
