from src.qodon.optimizers.optimizer import Optimizer
import numpy as np
import tensorflow as tf
import tensorflow_probability as tfp
from typing import List

class TfDiffEv(Optimizer):
    '''
    Tensorflow Differential Evolution optimizer for codon optimization.

    '''
    def __init__(self, config):
        super().__init__(config)
        self._optimize()

    def _optimize(self):
        '''
        Main execution. Run tensorflow optimizer. Objective
        function computes RNA structure with D-Wave's SA algorithm.

        '''

        self.initial_members = tf.convert_to_tensor(([_ for _ in self.config.initial_sequences]),np.float32)

        # Differential_weight: controls strength of mutations. We basically want to turn this off.
        # Crossover_prob: set this low. Need to think more about why this helps.
        optim_results = tfp.optimizer.differential_evolution_minimize(
            self._objective,
            initial_population=self.initial_members,
            max_iterations=self.config.args.codon_iterations,
            differential_weight=0.01,
            crossover_prob=0.1,
        )

        # Assign results as class attributes
        self.final_population = self._convert_to_ints(optim_results.final_population)
        self.final_energies = optim_results.final_objective_values
        self.mfe = np.min(optim_results.final_objective_values)
        self.mfe_index = np.argmin(optim_results.final_objective_values)

        self.final_codons = self._reverse_translate(self.final_population)
        self._verify_dna(self.final_codons[self.mfe_index])

        print(self.mfe)
        print(self.final_codons[self.mfe_index])

    def _objective(self, members):
        '''
        Objective function for TF to minimize

        NOTE: TF uses gradient descent to minimize continuous valued functions.
        The approach used here is not mathematically sound. It's a hack. But
        it gets the job done.

        '''

        # Map continuous valued tensor to RNA sequence
        n_seqs = self._convert_to_ints(members)
        n_seqs = self._reverse_translate(n_seqs)

        # Use the imported scoring function to score all sequences.
        scores = [self._tf_fold(s) for s in n_seqs]

        # Return TF object
        return tf.cast(scores, np.float32)

    def _convert_to_ints(self, members) -> List:
        '''
        Continuous --> discrete transformation

        Doesn't make mathematical sense but it works.

        '''
        # This is a hack. TF deals with continuous valued functions. We need discrete and finite.
        # So let's cheat. Whatever values are assigned, make them ints and take the absolute value.
        members = np.absolute(np.array(members).astype(int))

        return members
