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
        self._post_process()

    def _optimize(self):
        """
        Main method for random codon optimization

        """
        if not self.config.args.resume:
            self._iterate(self.initial_sequences, False)

        for i in range(self.config.args.codon_iterations):
            extra_sequences = self._generate_sequences(self.config.args.n_trials)
            self._iterate(extra_sequences)
