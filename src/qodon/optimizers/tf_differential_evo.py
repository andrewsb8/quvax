from src.qodon.optimizer import CodonOptimizer
import numpy as np
from typing import List
import tensorflow as tf
import tensorflow_probability as tfp


class TfDiffEv(CodonOptimizer):
    """
    Tensorflow Differential Evolution optimizer for codon optimization.

    """

    def __init__(self, config):
        super().__init__(config)
        tf.random.set_seed(self.config.args.random_seed)
        # tensorflow counts initial pop as first step but others don't
        # line below makes counting consistent among all optimizers
        self.codon_optimize_step -= 1
        #initialize tensor to store energies and list to store its indices
        #see Scalar Updates here: https://www.tensorflow.org/api_docs/python/tf/tensor_scatter_nd_update
        self.energies_tensor = tf.Variable([0 for i in range(self.config.args.n_trials)], dtype=np.float32)
        self._optimize()
        self._post_process()

    def _optimize(self):
        """
        Main execution. Run tensorflow optimizer. Objective
        function computes RNA structure with D-Wave's SA algorithm.

        """

        if self.config.args.resume:
            members = [self._convert_codons_to_ints(s) for s in self.initial_sequences]
            members = tf.convert_to_tensor(([_ for _ in members]), np.float32)
        else:
            members = tf.convert_to_tensor(
                ([_ for _ in self.initial_sequences]), np.float32
            )

        # Differential_weight: controls strength of mutations. We basically want to turn this off.
        # Crossover_prob: set this low. Need to think more about why this helps.
        tfp.optimizer.differential_evolution_minimize(
            self._objective,
            initial_population=members,
            max_iterations=self.config.args.codon_iterations,
            differential_weight=0.01,
            crossover_prob=0.1,
            func_tolerance=-1,  # force tensorflow to do max_iterations
        )

    def _objective(self, members):
        """
        Objective function for TF to minimize

        NOTE: TF uses gradient descent to minimize continuous valued functions.
        The approach used here is not mathematically sound. It's a hack. But
        it gets the job done.

        """

        # Map continuous valued tensor to integers associated with RNA sequence
        n_seqs = self._convert_to_ints(members)
        self._iterate(n_seqs)

        # Return TF object
        return self.energies_tensor.assign(self.energies)

    def _convert_to_ints(self, members) -> List:
        """
        Continuous --> discrete transformation

        Doesn't make mathematical sense but it works.

        """
        # This is a hack. TF deals with continuous valued functions. We need discrete and finite.
        # So let's cheat. Whatever values are assigned, make them ints and take the absolute value.
        members = np.absolute(np.array(members).astype(int))

        return members
