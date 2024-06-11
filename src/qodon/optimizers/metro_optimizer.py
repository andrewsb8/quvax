from src.qodon.optimizer import CodonOptimizer
import math
import random
import numpy as np

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
        #number of proposed changes
        num_changes = 1
        #beta = 1/kT
        beta = 1
        energies = self.energies
        sec_structs = self.sec_structs

        for i in range(self.config.args.codon_iterations):
            for j, sequence in enumerate(members):
                #first propose a change in codon with our perturb function
                proposed_members = self._perturb_dna(sequence, num_changes)
                self._fold_rna(self._convert_ints_to_codons(proposed_members))
                #If the new energy is lower than the old energy, we will accept the proposed sequence.
                if self.folder.best_score <= energies[j]:
                    members[j] = proposed_members
                    energies[j] = self.folder.best_score
                    sec_structs[j] = self.folder.dot_bracket
                    accepted += 1
                #Otherwise, we need to generate a probability
                elif math.e**(-beta*(self.folder.best_score - energies[j])) >= random.uniform(0.0, 1.0):
                    members[j] = proposed_members
                    energies[j] = self.folder.best_score
                    sec_structs[j] = self.folder.dot_bracket
                    accepted += 1
                #Rejects the change if not, no reassignment necessary
                else:
                    rejected += 1
            self._update_codon_step()
            self._iterate(members, energies, sec_structs) #pass energies and ss to _iterate to avoid refolding
        self.config.log.info("Amount accepted: " + str(accepted))
        self.config.log.info("Amount rejected: " + str(rejected))

    def _perturb_dna(self, old_genes: list, num_changes):
        """
        Function that randomly changes codons (number of changes is num_changes)
        in the input sequence. Inputs the amino acid sequence, the codon indexes
        (old_genes), and the number of changes to be made. Returns changed codon
        indexes (same format as initial_sequences).

        """
        #first, randomly select num_change number of indices in the sequence
        indices = random.sample(range(0, len(self.config.seq)-1), num_changes)

        for k in range(len(indices)):
            #Define what amino acid is being changed
            change_res = self.config.seq[indices[k]]
            #Determine the old codon number so that we can ensure it is actually changed
            old_codon = old_genes[indices[k]]
            new_codon = old_codon
            num_codons = len(self.code_map[change_res]["codons"])
            #Use the code map to randomly change the codon
            while new_codon == old_codon:
                old_genes[indices[k]] = random.randint(0,num_codons)
                new_codon = old_genes[indices[k]]
        return old_genes
