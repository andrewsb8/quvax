from src.qodon.optimizer import CodonOptimizer


class RandomOptimizer(CodonOptimizer):
    """
    A basic implementation of a random optimizer for a codon sequence.

    Parameters
    ----------

    """

    def __init__(self, config):
        super().__init__(config)

    def _optimize(self):
        """
        Main method for random codon optimization

        """
        if not self.config.args.resume:
            self._iterate(self.initial_sequences, update_counter=False)

        for i in range(self.config.args.codon_iterations):
            self.members = self._generate_sequences(self.config.args.population_size)
            self._iterate()
        self._post_process()
