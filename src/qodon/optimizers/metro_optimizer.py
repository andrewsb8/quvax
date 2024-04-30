from src.qodon.optimizer import CodonOptimizer
import math
import random
import numpy as np
from copy import deepcopy

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
 
        """
        
    
        """


        #initializing energy array
        energy_array = np.zeros(self.config.args.n_trials)
        #initializing the population
        members = self.initial_sequences
        #initializing accepted and rejected counters
        accepted = 0
        rejected = 0
           
             
        for i in range(self.config.args.codon_iterations):
            #generating the random number of changes we will make at this step                                                 
            num_changes = random.randint(1,int(len(members[0])//4))
            print("Number of changes for codon step ",i, ": ", num_changes)
            for j, sequence in enumerate(members):    
                #first we want to propose a change in codon with our perturb function
                proposed_members = members[j].copy()
                proposed_members = self._perturb_dna(proposed_members, num_changes)    
                print("current indexes: ", members[j])
                print("proposed indexes: ", proposed_members)
                #Now, we want to get the energy of the fold for both sequences
                current_E = self._fold_rna(self._reverse_translate(members[j]))
                proposed_E = self._fold_rna(self._reverse_translate(proposed_members))
                #If the new energy is lower than the old energy, we will accept the proposed sequence.
                if proposed_E <= current_E:
                    members[j] = proposed_members
                    energy_array[j] = proposed_E
                    accepted += 1
                #Otherwise, we need to generate a probability
                else: 
                    dE = proposed_E - current_E #energy difference. This will always be positive.
                    T = 200
                    kb = .1 #The values of T and kb are arbitrary because the folding energy is not in units of joules
                    Beta = 1/(kb*T)
                    Prob = math.e**(-Beta*dE)
                    number = random.random()
                    #Accept the change at the rate given by the probability
                    if Prob >= number:
                        members[j] = proposed_members
                        energy_array[j] = proposed_E
                        accepted += 1
                    #Rejects the change at 1-Prob
                    else:
                        members[j] = members[j]
                        energy_array[j] = current_E
                        rejected += 1
            print("Iteration number, :", i, " Energy Array: ", energy_array)
        print("Amount accepted: ", accepted)
        print("Amount rejected: ", rejected) 
        return energy_array #this is an array that shows how the energy progresses

    def _perturb_dna(self, old_genes: list, num_changes):
        """

        Function that randomly changes one codon in the input sequence. 
        
        Inputs the amino acid sequence, the codon indexes (old_genes), and the number of changes to be made.

        Returns changed codon indexes (same format as initial_sequences).

        """
        #first, randomly select an index in the sequence
        sequence = self.config.seq #this is a string of amino acids
        indeces = []
        i = 0
        while i < num_changes:
            random_index = random.randint(0,len(sequence)-1)
            if random_index in indeces:
                i+=0
            else:
                indeces.append(random_index)
                i+=1
        for k in range(0,len(indeces)):
            #Define what amino acid is being changed
            change_res = sequence[indeces[k]]
            #Determine the old codon number so that we can ensure it is actually changed
            old_codon = old_genes[indeces[k]]
            new_codon = deepcopy(old_codon)
        #Use the code map to randomly change the codon
            while new_codon == old_codon:
                num_codons = len(self.code_map[change_res]["codons"])
                old_genes[indeces[k]] = random.randint(0,num_codons)
                new_codon = old_genes[indeces[k]]
                
        return old_genes

        
    
