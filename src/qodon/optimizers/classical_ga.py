from src.qodon.optimizer import CodonOptimizer
import random


class GeneticAlgorithm(CodonOptimizer):
    """
    A basic implementation of a genetic algorithm.

    Parameters
    ----------

    """

    def __init__(self, config):
        super().__init__(config)
        self._optimize()
        self._post_process()

    def __repr__(self):
        return "Classical genetic algorithm for codon optimization."

    def _optimize(self):
        """
        Main method for codon optimization

        """
        if not self.config.args.resume:
            self._iterate(self.initial_sequences, False)
            members = self.initial_sequences
        else:
            members = [self._convert_codons_to_ints(s) for s in self.initial_sequences]

        # Simulate evolution for number of codon_iterations specified by user
        for i in range(self.config.args.codon_iterations):
            # Introduce mutations
            members = self._procreate(members)
            self._iterate(members)

    def _procreate(self, eligible_members):
        """
        Simulate procreation by randomly picking two genes
        and randomly recombining them with mutations.

        """

        new_members = []
        for i_trial in range(len(eligible_members)):
            lucky_pair = random.sample(eligible_members, 2)
            new_members.append(
                self._mutate_dna(
                    self._mix_genes(lucky_pair[0], lucky_pair[1]), mutation_chance=0.05
                )
            )
        return new_members

    def _mix_genes(self, genes_xx, genes_xy):
        """
        Create new genes by randomly mixing two

        """
        new_genes = []
        for i in range(len(genes_xx)):
            random_chance = random.uniform(0.0, 1.0)
            if random_chance < 0.5:
                new_genes.append(genes_xx[i])
            else:
                new_genes.append(genes_xy[i])
        return new_genes

    # mutate_chance should be defined elsewhere, probably yaml
    def _mutate_dna(self, old_genes: list, mutation_chance=0.01):
        """
        Randomly introduce mutations

        """
        new_d_sequence = ""
        new_indices = []
        total_log_score = 0.0
        for i, res in enumerate(self.config.seq):
            if mutation_chance > random.uniform(0.0, 1.0):
                passing_indices = []
                for j, chance in enumerate(self.code_map[res]["probs"]):
                    if chance > random.uniform(0.0, 1.0):
                        passing_indices.append(j)
                chosen_index = passing_indices[0]
            else:
                chosen_index = old_genes[i]
            new_indices.append(chosen_index)
            total_log_score += self.code_map[res]["log_scores"][chosen_index]
            new_d_sequence += self.code_map[res]["codons"][chosen_index]
        return new_indices
