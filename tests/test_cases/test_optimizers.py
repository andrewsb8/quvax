import sys
import os
import unittest
from unittest.mock import patch
from src.params.parser import Parser

class TestOptimizers(unittest.TestCase):
    """
    Test codon optimization will run for each top level optimizer option

    """
    def test_TFDE(self):
        from src.qodon.optimizers.tf_differential_evo import TfDiffEv
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-n", "4", "-c", "4", "-ms", "2", "-co", "TFDE"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = TfDiffEv(parser)
            with self.assertLogs(level='INFO'):
                opt._verify_dna(opt.final_codons[opt.mfe_index])


    def test_RAND(self):
        from src.qodon.optimizers.random_optimizer import RandomOptimizer
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-n", "4", "-c", "4", "-ms", "2", "-co", "RAND"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = RandomOptimizer(parser)
            with self.assertLogs(level='INFO'):
                opt._verify_dna(opt.final_codons[opt.mfe_index])

    def test_GA(self):
        from src.qodon.optimizers.classical_ga import GeneticAlgorithm
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-n", "4", "-c", "4", "-ms", "2", "-co", "GA"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = GeneticAlgorithm(parser)
            with self.assertLogs(level='INFO'):
                opt._verify_dna(opt.final_codons[opt.mfe_index])

if __name__ == '__main__':
    unittest.main()
