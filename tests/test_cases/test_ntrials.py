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
            parser.config.codon_table, parser.config.codon_scores, parser.config.code_map = CodonOptimizer._construct_codon_table(parser, parser.args.species)
            initial_sequences = CodonOptimizer._get_initial_sequences(parser, parser.seq, parser.args.n_trials)
            self.assertEqual(len(initial_sequences), 5)

    def test_lenNTrials_str(self):
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-n", "2"]
        with patch.object(sys, 'argv', testargs):
            with self.assertRaises(TypeError):
                parser = Parser()
                parser.args.n_trials = str(parser.args.n_trials)
                self.codons = CodonOptimizer._get_initial_sequences(parser, parser.seq, parser.args.n_trials)

if __name__ == '__main__':
    unittest.main()
