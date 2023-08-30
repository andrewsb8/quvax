from qodon.codon_tables import code_map
import random

class GenerateInitialSequences(object):
    def __init__(self, seq, ntrials):
        self.ntrials = ntrials
        self.code_map = code_map
        self.initial_sequences = self._get_initial_sequences(seq)

    def _get_initial_sequences(self, seq):
        '''
        Initialize population with randomly assembled members.

        '''
        initial_members = []
        for i in range(self.ntrials):
            d_sequence = ""
            chosen_indices = []
            for res in seq:
                random_prob = random.uniform(0.0, 1.0)
                reference_chances = code_map[res]['probs']
                passing_indices = []
                for chance in reference_chances:
                    if chance > random_prob:
                        passing_indices.append(reference_chances.index(chance))
                chosen_index = passing_indices[0]
                chosen_indices.append(chosen_index)
                d_sequence += code_map[res]['codons'][chosen_index]
            #0 was the codon sequence score. Breaks if value is removed. Not sure why yet
            member = [0, chosen_indices]
            initial_members.append(member)
        return initial_members
