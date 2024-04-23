from src.qodon.optimizer import CodonOptimizer
import math
import random

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
        
        n_seqs = [self._reverse_translate(s) for s in self.initial_sequences] #storing the nucleotide sequences
        energies = [self._fold_rna(s) for s in n_seqs] #getting the folding energies of these sequences
        self._write_output(n_seqs, energies, None)
 
        print("self.config.seq: ", self.config.seq)
        print("initial_sequences: ", self.initial_sequences)
        print("n_seqs: ", n_seqs)
        """
        
        All of this is taken from the random optimizer script so far.

        """
        for j, s in enumerate(self.initial_sequences):
            current_sequence = s
            energy_array = []
            for i in range(self.config.args.codon_iterations):
                #first we want to propose a change in codon with our perturb function
                proposed_sequence = _perturb_dna(self, current_sequence)
                current_nucleo = self._reverse_translate(current_sequence)
                proposed_nucleo = self._reverse_translate(proposed_sequence)
                #Now, we want to get the energy of the fold for both sequences
                current_E = _fold_rna(self, current_nucleo)
                proposed_E = _fold_rna(self, proposed_nucleo)
                #If the new energy is lower than the old energy, we will accept the proposed sequence.
                if proposed_E < current_E:
                    current_sequence = proposed_sequence
                    energy_array[j,i] = proposed_E
                #Otherwise, we need to generate a probability
                else: 
                    dE = proposed_E - current_E #energy difference. This will always be positive.
                    T = 2
                    kb = 1 #need to figure out a good value for this as well as Temperature.
                    Beta = 1/(kb*T)
                    Prob = math.e**(-Beta*dE)
                    number = random.random()
                    #Accept the change at the rate given by the probability
                    if Prob >= number:
                        current_sequence = proposed_sequence
                        energy_array[j,i] = proposed_E
                    #Rejects the change at 1-Prob
                    else:
                        current_sequence = current_sequence
                        energy_array[j,i] = current_E
                        
        return energy_array #this is an array that shows how the energy progresses

    def _perturb_dna(self, old_genes: list):
        """

        Function that randomly changes one codon in the input sequence. 
        
        Inputs the amino acid sequence and the codon indexes (old_genes).

        Returns changed codon indexes.

        """
        #first, randomly select an index in the sequence
        sequence = self.config.seq #this is a string of amino acids
        random_index = random.randint(0,len(sequence))
        #Define what amino acid is being changed
        change_res = sequence[random_index]
        #Determine the old codon number so that we can ensure it is actually changed
        old_codon = old_genes[random_index]
        new_codon = old_codon
        #Use the code map to randomly change the codon
        i = 1
        while new_codon == old_codon:
            num_codons = len(self.code_map[change_res]["codons"])
            old_genes[random_index] = random.randint(0,num_codons)
            new_codon = old_genes[random_index]
        new_genes = old_genes
        return new_genes
        
        
