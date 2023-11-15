import sys
import os
import unittest
from unittest.mock import patch
from src.params.parser import Parser
from src.qodon.optimizer import CodonOptimizer
from Bio import SeqIO

class FakeOptimizer(CodonOptimizer):
    def __init__(self, config):
        super().__init__(config)

    def _optimize(self):
        pass

class TestNTrials(unittest.TestCase):
    """
    Test variable n_trials will generate appropriate number of sequences

    """
    def test_lenNTrials_5(self):
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-n", "5"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = FakeOptimizer(parser)
            self.assertEqual(len(opt.initial_sequences), 5)

    def test_lenNTrials_str(self):
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-n", "2"]
        with patch.object(sys, 'argv', testargs):
            with self.assertRaises(TypeError):
                parser = Parser()
                parser.args.n_trials = str(parser.args.n_trials)
                opt = FakeOptimizer(parser)

if __name__ == '__main__':
    unittest.main()
