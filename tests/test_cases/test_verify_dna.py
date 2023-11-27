import sys
import os
import unittest
from unittest.mock import patch
from src.params.parser import Parser
from tests.mock_optimizer import MockOptimizer

class TestVerifyDNA(unittest.TestCase):
    """
    Test Optimizer Method _verify_dna for successful and nonsuccessful translations

    """
    def test_NoTranslate(self):
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-n", "4", "-c", "4", "-ms", "2", "-co", "TFDE"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = MockOptimizer(parser)
            opt.final_codons = "ATGTTCGTATTCTTAGTGTTACTGCCGCTCGT"
            with self.assertRaises(ValueError):
                opt._verify_dna(opt.final_codons)

    def test_WillTranslate(self):
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-n", "4", "-c", "4", "-ms", "2", "-co", "TFDE"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = MockOptimizer(parser)
            opt.final_codons = "GGCGGCGGGAAC"
            with self.assertLogs(level='INFO'):
                opt._verify_dna(opt.final_codons)

if __name__ == '__main__':
    unittest.main()
