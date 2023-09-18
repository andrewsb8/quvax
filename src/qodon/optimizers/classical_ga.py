from src.qodon.optimizers.optimizer import Optimizer
from qodon.src.codon_tables import code_map
import random
from operator import itemgetter
from Bio.Seq import Seq


class GeneticAlgorithm(Optimizer):
    """
    A basic implementation of a genetic algorithm.

    Parameters
    ----------

    """
    def __init__(self, seq):
        self.seq = seq
        self.code_map = code_map
        self.elitelist = 10
        self.randomlist = 2
        self.ntrials = 10
        self.numgens = 100
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
        population = self._get_initial_members()

        # Simulate evolution for 'numgens' trials
        for i in range(self.numgens):

            # Sort sequences in ascending order by score
            ranked_members = sorted(population, key=itemgetter(0))

            # Isolate subset of sequences with best score
            fittest_members = ranked_members[:self.elitelist]

            # Randomly sample the remaining members
            lucky_members = random.sample(ranked_members[self.elitelist:],
                                                        self.randomlist)

            # Members eligible for mutation are 'best' and 'lucky'
            eligible_members = fittest_members + lucky_members

            # Introduce mutations
            population += self._procreate(eligible_members)

        # Record fittest member of population after simulating evo
        fittest_member = sorted(population, key=itemgetter(0))[0]
        self.optimal_codon_indices = fittest_member[1]

    def _procreate(self, eligible_members):
        '''
        Simulate procreation by randomly picking two genes
        and randomly recombining them with mutations.

        '''

        new_members = []
        for i_trial in range(self.ntrials):
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

    def _mutate_dna(self, old_genes: list, mutation_chance=0.01):
        '''
        Randomly introduce mutations

        '''
        new_d_sequence = ""
        new_indices = []
        total_log_score = 0.0
        for i, res in enumerate(self.seq):
            if mutation_chance > random.uniform(0.0, 1.0):
                passing_indices = []
                for j, chance in enumerate(self.code_map[res]['probs']):
                    if chance > random.uniform(0.0, 1.0):
                        passing_indices.append(j)
                chosen_index = passing_indices[0]
            else:
                chosen_index = old_genes[i]
            new_indices.append(chosen_index)
            total_log_score += self.code_map[res]['log_scores'][chosen_index]
            new_d_sequence += self.code_map[res]['codons'][chosen_index]
        #0 was the codon sequence score. Breaks if value is removed. Not sure why yet
        return [0, new_indices]

    def _reverse_translate(self):
        '''
        Convert to nucleotide sequence

        '''
        self.n_seq = ''.join([
            self.code_map[res]['codons'][self.optimal_codon_indices[i]]
            for i, res in enumerate(self.seq)
        ])

    def _verify_dna(self):
        '''
        Translate nucleotide sequence to make sure it matches input

        '''
        if self.seq != str(Seq(self.n_seq).transcribe().translate()):
            raise ValueError(
                "Error: Codon sequence did not translate properly!")
