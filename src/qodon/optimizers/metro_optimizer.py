from src.qodon.optimizer import CodonOptimizer

class MetropolisOptimizer(CodonOptimizer):

    """
    Implementation of the Monte Carlo Metropolis Algorithm for a codon sequence

    Parameters

    -----------

    """
    
    def __init__(self, config):
        super().__init__(config)
        self.optimize()

    def _optimize(self):
        """

        Method for codon optimization using metropolis algorithm
        
        """
        
        n_seqs = [self.reverse_translate(s) for s in self.initial_sequences] #storing the nucleotide sequences
        energies = [self._fold_rna(s) for s in n_seqs] #getting the folding energies of these sequences
        self._extend_output(n_seqs, energies, None)
        
        """
        
        All of this is taken from the random optimizer script so far.

        """
    
        for i in range(self.config.args.codon_iterations):
            
