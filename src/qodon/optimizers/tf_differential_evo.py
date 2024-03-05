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
        self._optimize()

    def _optimize(self):
        """
        Main execution. Run tensorflow optimizer. Objective
        function computes RNA structure with D-Wave's SA algorithm.

        """

        self.initial_members = tf.convert_to_tensor(
            ([_ for _ in self.initial_sequences]), np.float32
        )

        # Differential_weight: controls strength of mutations. We basically want to turn this off.
        # Crossover_prob: set this low. Need to think more about why this helps.
        tfp.optimizer.differential_evolution_minimize(
            self._objective,
            initial_population=self.initial_members,
            max_iterations=self.config.args.codon_iterations,
            differential_weight=0.01,
            crossover_prob=0.1,
            func_tolerance=-1, #force tensorflow to do max_iterations
        )

        self._get_optimized_sequences()
        if self.config.args.target is not None:
            self._check_target()
        self._pickle_output()

    def _objective(self, members):
        """
        Objective function for TF to minimize

        NOTE: TF uses gradient descent to minimize continuous valued functions.
        The approach used here is not mathematically sound. It's a hack. But
        it gets the job done.

        """

        # Map continuous valued tensor to RNA sequence
        n_seqs = self._convert_to_ints(members)
        n_seqs = [self._reverse_translate(s) for s in n_seqs]

        # Use the imported scoring function to score all sequences.
        energies = [self._fold_rna(s) for s in n_seqs]

        self._extend_output(n_seqs, energies, None)

        # Return TF object
        return tf.cast(energies, np.float32)

    def _convert_to_ints(self, members) -> List:
        """
        Continuous --> discrete transformation

        Doesn't make mathematical sense but it works.

        """
        # This is a hack. TF deals with continuous valued functions. We need discrete and finite.
        # So let's cheat. Whatever values are assigned, make them ints and take the absolute value.
        members = np.absolute(np.array(members).astype(int))

        return members
