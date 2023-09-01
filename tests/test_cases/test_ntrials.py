import sys
sys.path.append('../../')
import os
import unittest
from qodon.initiate_sequences import GenerateInitialSequences

class TestNTrials(unittest.TestCase):
    def test_lenNTrials_5(self):
        self.seq = "GGG"
        self.ntrials = 5
        self.codons = GenerateInitialSequences(self.seq, self.ntrials)
        self.assertEqual(len(self.codons.initial_sequences), 5)

    def test_lenNTrials_10(self):
        self.seq = "GGG"
        self.ntrials = 10
        self.codons = GenerateInitialSequences(self.seq, self.ntrials)
        self.assertEqual(len(self.codons.initial_sequences), 10)

if __name__ == '__main__':
    unittest.main()
