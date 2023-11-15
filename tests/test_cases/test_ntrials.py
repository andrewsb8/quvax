import sys
import os
import unittest
from unittest.mock import patch
from src.params.parser import Parser
from src.qodon.optimizer import CodonOptimizer
from Bio import SeqIO

class TestNTrials(unittest.TestCase):
    """
    Test variable n_trials will generate appropriate number of sequences

    """
    def test_lenNTrials_5(self):
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-n", "5"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = CodonOptimizer(parser)
            opt.codon_table, opt.codon_scores, opt.code_map = opt._construct_codon_table(opt)
            opt.initial_sequences = opt._get_initial_sequences(opt)
            self.assertEqual(len(opt.initial_sequences), 5)

    def test_lenNTrials_str(self):
        self.seq = str(SeqIO.read('tests/test_sequences/GGGN.fasta','fasta').seq)
        self.n_trials = "2"
        with self.assertRaises(TypeError):
            self.codons = GenerateInitialSequences(self.seq, self.n_trials)

if __name__ == '__main__':
    unittest.main()
