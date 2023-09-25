from src.qodon.optimizers.optimizer import Optimizer
from src.qodon.codon_tables import code_map
import random
from operator import itemgetter
import numpy as np


class GeneticAlgorithm(Optimizer):
    """
    A basic implementation of a genetic algorithm.

    Parameters
    ----------

    """
    def __init__(self, config):
        super().__init__(config)

        #parameters are used in _propagate_generations() - should they be added to yaml parameters? removed?
        self.elitelist = 10
        self.randomlist = 2

        self._optimize()

    def __repr__(self):
        return 'Classical genetic algorithm for codon optimization.'

    def _optimize(self):
        '''
        Main method for codon optimization

        '''

        # Simulate evolution
        self._propagate_generations()

        # Recover nucleotide sequence
        self._reverse_translate()

        # Make sure fittest member translates to correct aa sequence
        self._verify_dna()

    def _propagate_generations(self):

        # Initialize population
        n_seqs = self.config.initial_sequences

        # Simulate evolution for number of codon_iterations specified by user
        for i in range(self.config.args.codon_iterations):
            # Introduce mutations
            n_seqs += self._procreate(n_seqs)
            # Use the imported scoring function to score all sequences.
            scores = [self._tf_fold(s) for s in n_seqs]

        # Record fittest member of population after simulating evo
        self.mfe = np.min(scores)
        self.final_codons = n_seqs[np.argmin(scores)]

    def _procreate(self, eligible_members):
        '''
        Simulate procreation by randomly picking two genes
        and randomly recombining them with mutations.

        '''

        new_members = []
        for i_trial in range(len(eligible_members)):
            lucky_pair = random.sample(eligible_members, 2)
            new_members.append(
                self._mutate_dna(self._mix_genes(lucky_pair[0][1],
                                                 lucky_pair[1][1]),
                                 mutation_chance=0.05))
        return new_members

    @staticmethod
    def _mix_genes(genes_xx, genes_xy):
        '''
        Create new genes by randomly mixing two

        '''
        new_genes = []
        for i in range(len(genes_xx)):
            random_chance = random.uniform(0.0, 1.0)
            if random_chance < 0.5:
                new_genes.append(genes_xx[i])
            else:
                new_genes.append(genes_xy[i])
        return new_genes

    #mutate_chance should be defined elsewhere, probably yaml
    def _mutate_dna(self, old_genes: list, mutation_chance=0.01):
        '''
        Randomly introduce mutations

        '''
        new_d_sequence = ""
        new_indices = []
        total_log_score = 0.0
        for i, res in enumerate(self.config.seq):
            if mutation_chance > random.uniform(0.0, 1.0):
                passing_indices = []
                for j, chance in enumerate(self.config.code_map[res]['probs']):
                    if chance > random.uniform(0.0, 1.0):
                        passing_indices.append(j)
                chosen_index = passing_indices[0]
            else:
                chosen_index = old_genes[i]
            new_indices.append(chosen_index)
            total_log_score += self.config.code_map[res]['log_scores'][chosen_index]
            new_d_sequence += self.config.code_map[res]['codons'][chosen_index]
        #0 was the codon sequence score. Breaks if value is removed. Not sure why yet
        return [0, new_indices]
