import sys
import os
import unittest
from src.qodon.initiate_sequences import GenerateInitialSequences
from Bio import SeqIO

class TestNTrials(unittest.TestCase):
    """
    Test variable ntrials will generate appropriate number of sequences

    """
    def test_lenNTrials_5(self):
        self.seq = str(SeqIO.read('tests/test_sequences/GGGN.fasta','fasta').seq)
        self.ntrials = 5
        self.codons = GenerateInitialSequences(self.seq, self.ntrials)
        self.assertEqual(len(self.codons.initial_sequences), 5)

    def test_lenNTrials_str(self):
        self.seq = str(SeqIO.read('tests/test_sequences/GGGN.fasta','fasta').seq)
        self.ntrials = "2"
        with self.assertRaises(TypeError):
            self.codons = GenerateInitialSequences(self.seq, self.ntrials)

if __name__ == '__main__':
    unittest.main()
