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
        n_seqs = [self._reverse_translate(s) for s in self.initial_sequences]
        energies = [self._fold_rna(s) for s in n_seqs]
        self._extend_output(n_seqs, energies, None)

        for i in range(self.config.args.codon_iterations):
            self._update_codon_step()
            extra_sequences = self._generate_sequences(self.config.args.n_trials)
            n_seqs = [self._reverse_translate(s) for s in extra_sequences]
            energies = [self._fold_rna(s) for s in n_seqs]
            self._extend_output(n_seqs, energies, None)
            
        self._get_optimized_sequences()
        if self.config.args.target is not None:
            self._check_target()
        self._pickle_output()
