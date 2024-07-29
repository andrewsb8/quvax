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
    minimum trapping, redundant sequences in the database, and increase sampling
    of phase space. The user-defined convergence criterion will still allow the
    optimization to terminate in a reasonable amount of iterations.

    Parameters
    -----------

    """

    def __init__(self, config):
        super().__init__(config)
        self.accepted = 0
        self.rejected = 0
        self.randomed = 0

    def _optimize(self):
        """
        Method for codon optimization using metropolis algorithm

        """

        for i in range(self.config.args.codon_iterations):
            self._metropolis_iteration(self.members, self.energies, self.sec_structs)
            self._iterate(fold_sequences=False) #do not calculated folding energy or secondary structures
        self.config.log.info(
            "MC stats are only kept track of in each individual execution of design.py and are not aggregated from previous runs if using --resume."
        )
        self.config.log.info("Sequence changes accepted: " + str(self.accepted))
        self.config.log.info("Sequence changes rejected: " + str(self.rejected))
        self.config.log.info("Random sequences generated: " + str(self.randomed))
        self._post_process()

    def _metropolis_iteration(self, members, energies, sec_structs):
        for j, sequence in enumerate(members):
            seq_copy = copy.deepcopy(sequence)
            self.seq_rejections = 0
            # try multiple changes in case of rejections
            while True:
                # if reach maximum rejects, randomly sample another sequence
                if self.seq_rejections >= self.config.args.sequence_rejections:
                    members[j] = self._generate_sequences(1)[0]
                    self._fold_rna(self._convert_ints_to_codons(members[j]))
                    energies[j] = self.folder.best_score
                    sec_structs[j] = self.folder.dot_bracket
                    self.randomed += 1
                    self.rejected += self.seq_rejections
                    break
                # propose change and fold
                proposed_members = self._perturb_dna(
                    seq_copy, self.config.args.num_sequence_changes
                )
                self._fold_rna(self._convert_ints_to_codons(proposed_members))
                # require sequence is not the same which is an edge case (i.e.
                # propose change to only amino acids with one possible codon)
                # that can save some cycles
                if proposed_members != sequence and self._check_changes(
                    self.folder.best_score,
                    self.config.args.beta,
                    self.energies[j],
                    self.config.args.beta,
                ):
                    self._accept_changes(
                        proposed_members, j, members, energies, sec_structs
                    )
                    break
                else:
                    self.seq_rejections += 1

    def _check_changes(self, energy, beta, ref_energy, ref_beta):
        # Accept lower energy,
        if beta * energy <= ref_beta * ref_energy:
            return True
        # Otherwise, we need to generate a probability
        elif self._calculate_boltzmann(
            energy, beta, ref_energy, ref_beta
        ) >= random.uniform(0.0, 1.0):
            return True
        # Rejects the change if not, no reassignment necessary
        else:
            return False

    def _accept_changes(self, proposed_members, index, members, energies, sec_structs):
        members[index] = proposed_members
        energies[index] = self.folder.best_score
        sec_structs[index] = self.folder.dot_bracket
        self.accepted += 1
        self.rejected += self.seq_rejections

    def _perturb_dna(self, old_genes: list, num_changes):
        """
        Function that randomly changes codons (number of changes is num_changes)
        in the input sequence. Inputs the amino acid sequence, the codon indexes
        (old_genes), and the number of changes to be made. Returns changed codon
        indexes (same format as initial_sequences).

        """
        # first, randomly select num_change number of indices in the sequence
        indices = random.sample(
            range(0, len(self.config.protein_sequence) - 1), num_changes
        )

        for k in range(len(indices)):
            # Define what amino acid is being changed
            change_res = self.config.protein_sequence[indices[k]]
            # Determine the old codon number so that we can ensure it is actually changed
            old_codon = old_genes[indices[k]]
            new_codon = old_codon
            num_codons = len(self.code_map[change_res]["codons"])
            if num_codons <= 1:
                return old_genes
            # Use the code map to randomly change the codon
            while new_codon == old_codon:
                # num_codons-1 because if randint returns num_codons, that is
                # equivalent to the 0th index and only num_changes-1 will be
                # truly applied, which could lead to errors
                old_genes[indices[k]] = random.randint(0, num_codons - 1)
                new_codon = old_genes[indices[k]]
        return old_genes

    def _calculate_boltzmann(self, energy, beta, ref_energy, ref_beta):
        return math.e ** (-((beta * energy) - (ref_beta * ref_energy)))
