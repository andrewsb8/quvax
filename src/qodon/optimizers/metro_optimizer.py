from src.qodon.optimizer import CodonOptimizer

class MetropolisOptimizer(CodonOptimizer):

    """
    Implementation of the Monte Carlo Metropolis Algorithm for a codon sequence

    Parameters

    -----------

    """
    
    def __init__(self, config):
        super().__init__(config)
        self._optimize()

    def _optimize(self):
        """

        Method for codon optimization using metropolis algorithm
        
        """
        
        n_seqs = [self.reverse_translate(s) for s in self.initial_sequences] #storing the nucleotide sequences
        energies = [self._fold_rna(s) for s in n_seqs] #getting the folding energies of these sequences
        self._extend_output(n_seqs, energies, None)
 
        print("n_seqs: ", n_seqs)
        print("initial_sequences: ", self.initial_sequences)

        """
        
        All of this is taken from the random optimizer script so far.

        """
    
#        for i in range(self.config.args.codon_iterations):
#            for s in n_seqs:
#                starting_sequence = n_seq[s]
                
                #now propose a change in codon sequence. Stuck here rn. What function/new code would I use to just change one codon? I am thinking possibly using _construct_codon_table, but not sure how to access the codon mappings. Might be some info in _generate_sequences that I can use.

                
