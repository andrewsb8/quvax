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
            self.n_seqs = [self._reverse_translate(s) for s in self.initial_sequences]
            self.energies = [self._fold_rna(s) for s in self.n_seqs]
            self._write_output(self.n_seqs, self.energies, None)

        for i in range(self.config.args.codon_iterations):
            extra_sequences = self._generate_sequences(self.config.args.n_trials)
            self.n_seqs = [self._reverse_translate(s) for s in extra_sequences]
            self.energies = [self._fold_rna(s) for s in self.n_seqs]
            self._update_mfe(self.energies)
            self._update_codon_step()
            self._write_output(self.n_seqs, self.energies, None)

        self._get_optimized_sequences()
        if self.config.args.target is not None:
            self._check_target()
