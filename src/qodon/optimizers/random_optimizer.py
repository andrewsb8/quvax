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
        Main method for codon optimization

        """

        num_extra_sequences = (
            self.config.args.codon_iterations * self.config.args.n_trials
        ) - self.config.args.n_trials
        extra_sequences = self._generate_sequences(num_extra_sequences)
        self.initial_sequences.extend(extra_sequences)
        n_seqs = [self._reverse_translate(s) for s in self.initial_sequences]

        energies = [self._fold_rna(s) for s in n_seqs]

        self._extend_output(n_seqs, energies, None)
        self._get_optimized_sequence()
        self._pickle_output()
