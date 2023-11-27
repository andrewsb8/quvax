import sys
import os
import unittest
from unittest.mock import patch
from src.params.parser import Parser
from src.qodon.optimizer import CodonOptimizer

class FakeOptimizer(CodonOptimizer):
    def __init__(self, config):
        super().__init__(config)

    def _optimize(self):
        pass

class TestVerifyDNA(unittest.TestCase):
    """
    Test codon optimization will run for each top level optimizer option

    """
    def test_NoTranslate(self):
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-n", "4", "-c", "4", "-ms", "2", "-co", "TFDE"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = FakeOptimizer(parser)
            opt.final_codons = "ATGTTCGTATTCTTAGTGTTACTGCCGCTCGT"
            with self.assertRaises(ValueError):
                opt._verify_dna(opt.final_codons)

    def test_WillTranslate(self):
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-n", "4", "-c", "4", "-ms", "2", "-co", "TFDE"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = FakeOptimizer(parser)
            opt.final_codons = "GGCGGCGGGAAC"
            with self.assertLogs(level='INFO'):
                opt._verify_dna(opt.final_codons)

if __name__ == '__main__':
    unittest.main()
