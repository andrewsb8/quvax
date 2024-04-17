from src.qodon.optimizer import CodonOptimizer


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
        """
        Main method for random codon optimization

        """
        if not self.config.args.resume:
            self._iterate(self.initial_sequences)

        for i in range(self.config.args.codon_iterations):
            extra_sequences = self._generate_sequences(self.config.args.n_trials)
            self._iterate(extra_sequences)
            self._update_codon_step()

        self._get_optimized_sequences()
        if self.config.args.target is not None:
            self._check_target()
