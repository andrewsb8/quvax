'''

    mRNA Codon Optimization with Quantum Computers
    Copyright (C) 2021  Dillion M. Fox, Ross C. Walker

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''
from qodon.src.codon_tables import code_map
import random

class GenerateInitialSequences(object):
    def __init__(self, seq):
        self.ntrials = 10
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
