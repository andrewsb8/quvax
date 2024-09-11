from src.codon_opt.codon_optimizer import CodonOptimizer
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
        # initialize tensor to store energies and list to store its indices
        # see Scalar Updates here: https://www.tensorflow.org/api_docs/python/tf/tensor_scatter_nd_update
        self.energies_tensor = tf.Variable(
            [0 for i in range(self.config.args.population_size)], dtype=np.float32
        )

    def _optimize(self):
        """
        Main execution. Run tensorflow optimizer. Objective
        function computes RNA structure with D-Wave's SA algorithm.

        """

        members = tf.convert_to_tensor(([_ for _ in self.members]), np.float32)

        tfp.optimizer.differential_evolution_minimize(
            self._objective,
            initial_population=members,
            max_iterations=self.config.args.codon_iterations,
            differential_weight=self.config.args.mutation_chance,
            crossover_prob=self.config.args.crossover_probability,
            func_tolerance=-1,  # force tensorflow to do max_iterations
        )
        self._post_process()

    def _objective(self, members):
        """
        Objective function for TF to minimize

        NOTE: TF uses gradient descent to minimize continuous valued functions.
        The approach used here is not mathematically sound. It's a hack. But
        it gets the job done.

        """

        # Map continuous valued tensor to integers associated with RNA sequence
        self.members = self._convert_to_ints(members)
        self._iterate()

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
