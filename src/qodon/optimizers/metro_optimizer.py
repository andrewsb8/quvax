from src.qodon.optimizer import CodonOptimizer
import math
import random
import numpy as np
import copy


class MetropolisOptimizer(CodonOptimizer):
    """
    Implementation of the Monte Carlo Metropolis Algorithm for a codon sequence.
    The algorithm will propose a user defined number of changes to each sequence
    in the population. If accepted, move to next sequence. If rejected, propose
    new changes. If a user-defined number of rejections happen for the same
    sequence, then a new random sequence is defined. This will prevent local
    minimum trapping and redundant sequences in the database. The user-defined
    convergence criterion will still allow the optimization to terminate in a
    reasonable amount of iterations.

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
            self._iterate(self.initial_sequences, update_counter=False)
            members = self.initial_sequences
        else:
            members = [self._convert_codons_to_ints(s) for s in self.initial_sequences]

        accepted = 0
        rejected = 0
        randomed = 0
        energies = self.energies
        sec_structs = self.sec_structs

        for i in range(self.config.args.codon_iterations):
            for j, sequence in enumerate(members):
                seq_copy = copy.deepcopy(sequence)
                seq_rejections = 0
                # try multiple changes in case of rejections
                while True:
                    # if reach maximum rejects, randomly sample another sequence
                    if seq_rejections >= self.config.args.sequence_rejections:
                        rand_seq = self._generate_sequences(1)
                        self._fold_rna(self._convert_ints_to_codons(rand_seq[0]))
                        energies[j] = self.folder.best_score
                        randomed += 1
                        rejected += seq_rejections
                        break
                    # propose change and fold
                    proposed_members = self._perturb_dna(
                        seq_copy, self.config.args.num_sequence_changes
                    )
                    self._fold_rna(self._convert_ints_to_codons(proposed_members))
                    # Accept lower energy, require sequence is not the same
                    # which is an edge case (i.e. propose change to only amino
                    # acids with one possible codon) that can save some cycles
                    if (
                        proposed_members != sequence
                        and self.folder.best_score <= energies[j]
                    ):
                        members[j] = proposed_members
                        energies[j] = self.folder.best_score
                        sec_structs[j] = self.folder.dot_bracket
                        accepted += 1
                        rejected += seq_rejections
                        break
                    # Otherwise, we need to generate a probability
                    elif proposed_members != sequence and math.e ** (
                        -self.config.args.beta * (self.folder.best_score - energies[j])
                    ) >= random.uniform(0.0, 1.0):
                        members[j] = proposed_members
                        energies[j] = self.folder.best_score
                        sec_structs[j] = self.folder.dot_bracket
                        accepted += 1
                        rejected += seq_rejections
                        break
                    # Rejects the change if not, no reassignment necessary
                    else:
                        seq_rejections += 1
            self._iterate(
                members, energies, sec_structs
            )  # pass energies and ss to _iterate to avoid refolding
        self.config.log.info("Sequence changes accepted: " + str(accepted))
        self.config.log.info("Sequence changes rejected: " + str(rejected))
        self.config.log.info("Random sequences generated: " + str(randomed))

    def _perturb_dna(self, old_genes: list, num_changes):
        """
        Function that randomly changes codons (number of changes is num_changes)
        in the input sequence. Inputs the amino acid sequence, the codon indexes
        (old_genes), and the number of changes to be made. Returns changed codon
        indexes (same format as initial_sequences).

        """
        # first, randomly select num_change number of indices in the sequence
        indices = random.sample(range(0, len(self.config.seq) - 1), num_changes)

        for k in range(len(indices)):
            # Define what amino acid is being changed
            change_res = self.config.seq[indices[k]]
            # Determine the old codon number so that we can ensure it is actually changed
            old_codon = old_genes[indices[k]]
            new_codon = old_codon
            num_codons = len(self.code_map[change_res]["codons"])
            if num_codons <= 1:
                return old_genes
            # Use the code map to randomly change the codon
            while new_codon == old_codon:
                old_genes[indices[k]] = random.randint(0, num_codons)
                new_codon = old_genes[indices[k]]
        return old_genes
