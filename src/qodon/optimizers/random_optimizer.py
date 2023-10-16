from src.qodon.optimizers.optimizer import Optimizer
from src.qodon.initiate_sequences import GenerateInitialSequences
from src.qodon.codon_tables import code_map
import random
from operator import itemgetter
import numpy as np


class RandomOptimizer(Optimizer):
    """
    A basic implementation of a genetic algorithm.

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

        extra_sequences = GenerateInitialSequences(self.config.seq, self.config.args.codon_iterations - self.config.args.n_trials).initial_sequences
        for sequence in extra_sequences:
            self.config.initial_sequences.append(sequence)

        scores = [self._tf_fold(s) for s in self.config.initial_sequences]

        self.final_energies = scores
        self.mfe = np.min(scores)
        self.mfe_index = np.argmin(scores)

        self.final_codons = self._reverse_translate(self.config.initial_sequences[self.mfe_index])

        self._verify_dna(self.final_codons)

        print(self.mfe)
        print(self.final_codons)
