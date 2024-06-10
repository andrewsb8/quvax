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
        self._post_process()

    def _optimize(self):
        """
        Method for codon optimization using metropolis algorithm

        """

        if not self.config.args.resume:
            self._iterate(self.initial_sequences)
            members = self.initial_sequences
        else:
            members = [self._convert_codons_to_ints(s) for s in self.initial_sequences]

        accepted = 0
        rejected = 0
        avg_changes = 0
        #number of proposed changes
        num_changes = 1
        #beta = 1/kT
        beta = 1

        for i in range(self.config.args.codon_iterations):
            for j, sequence in enumerate(members):
                #first propose a change in codon with our perturb function
                proposed_members = self._perturb_dna(members[j], num_changes)
                #Get the energy of the fold for both sequences
                current_E = self._fold_rna(self._convert_ints_to_codons(members[j]))
                proposed_E = self._fold_rna(self._convert_ints_to_codons(proposed_members))
                #If the new energy is lower than the old energy, we will accept the proposed sequence.
                if proposed_E <= current_E:
                    members[j] = proposed_members
                    energies[j] = proposed_E
                    accepted += 1
                #Otherwise, we need to generate a probability
                elif math.e**(-beta*(proposed_E - current_E)) >= random.uniform(0.0, 1.0):
                    members[j] = proposed_members
                    energies[j] = proposed_E
                    accepted += 1
                #Rejects the change if not
                else:
                    members[j] = members[j]
                    energies[j] = current_E
                    rejected += 1
            self._update_codon_step()
            self._iterate(members)
        print("Amount accepted: ", accepted)
        print("Amount rejected: ", rejected)

    def _perturb_dna(self, old_genes: list, num_changes):
        """
        Function that randomly changes one codon in the input sequence. Inputs
        the amino acid sequence, the codon indexes (old_genes), and the number
        of changes to be made. Returns changed codon indexes (same format
        as initial_sequences).

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
