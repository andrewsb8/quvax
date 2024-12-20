from src.codon_opt.codon_optimizer import CodonOptimizer
import random


class GeneticAlgorithm(CodonOptimizer):
    """
    A basic implementation of a genetic algorithm for codon optimization.

    Parameters
    ----------

    """

    def __init__(self, config):
        super().__init__(config)

    def _optimize(self):
        """
        Main method for codon optimization

        """

        # Simulate evolution for number of codon_iterations specified by user
        for i in range(self.config.args.codon_iterations):
            # Introduce mutations
            self.members = self._procreate(self.members)
            self._iterate()
        self._post_process()

    def _procreate(self, eligible_members):
        """
        Simulate procreation by randomly picking two genes
        and randomly recombining them with mutations.

        """

        new_members = []
        for i_trial in range(len(eligible_members)):
            lucky_pair = random.sample(eligible_members, 2)
            new_members.append(
                self._mutate_dna(self._mix_genes(lucky_pair[0], lucky_pair[1]))
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

    def _mutate_dna(self, old_genes: list):
        """
        Randomly introduce mutations

        """
        new_d_sequence = ""
        new_indices = []
        total_log_score = 0.0
        for i, res in enumerate(self.config.protein_sequence):
            if self.config.args.mutation_chance > random.uniform(0.0, 1.0):
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
