import sys
import os
import unittest
from unittest.mock import patch
from src.params.parser import Parser

class TestRandomSeed(unittest.TestCase):
    """
    Test random seed for initial sequence generation and optimization processes

    """
    def test_SameInitialSequences_DefaultSeed(self, _mock_optimizer):
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-n", "4", "-sd", "1"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = _mock_optimizer(parser)
            opt2 = _mock_optimizer(parser)
            self.assertEqual(opt.initial_sequences, opt2.initial_sequences)

    def test_SameInitialSequences_NewSeed(self, _mock_optimizer):
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-n", "4", "-sd", "1098760354"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = _mock_optimizer(parser)
            opt2 = _mock_optimizer(parser)
            self.assertEqual(opt.initial_sequences, opt2.initial_sequences)

    def test_TfDiffEv_DefaultSeed(self):
        from src.qodon.optimizers.tf_differential_evo import TfDiffEv
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-n", "4", "-c", "1000", "-sd", "1"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = TfDiffEv(parser)
            opt2 = TfDiffEv(parser)
            self.assertEqual(opt.final_codons, opt2.final_codons)

    def test_TfDiffEv_NewSeed(self):
        from src.qodon.optimizers.tf_differential_evo import TfDiffEv
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-co", "TFDE", "-n", "4", "-c", "1000", "-sd", "2546345746583"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = TfDiffEv(parser)
            opt2 = TfDiffEv(parser)
            self.assertEqual(opt.final_codons, opt2.final_codons)

    def test_TfDiffEv_OneIteration(self):
        from src.qodon.optimizers.tf_differential_evo import TfDiffEv
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-co", "TFDE", "-n", "4", "-c", "1", "-sd", "1"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = TfDiffEv(parser)
            opt2 = TfDiffEv(parser)
            self.assertEqual(opt.optimization_process['sequences'][-1], opt2.optimization_process['sequences'][-1])

    def test_TfDiffEv_LongerSequence(self):
        from src.qodon.optimizers.tf_differential_evo import TfDiffEv
        testargs = ["design.py", "-i", "examples/spike_trim_20.fasta", "-co", "TFDE", "-n", "4", "-c", "1000", "-sd", "1"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = TfDiffEv(parser)
            opt2 = TfDiffEv(parser)
            self.assertEqual(opt.final_codons, opt2.final_codons)

    def test_GA_DefaultSeed(self):
        from src.qodon.optimizers.classical_ga import GeneticAlgorithm
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-co", "GA", "-n", "4", "-c", "1000", "-sd", "1"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = GeneticAlgorithm(parser)
            opt2 = GeneticAlgorithm(parser)
            self.assertEqual(opt.final_codons, opt2.final_codons)

    def test_GA_NewSeed(self):
        from src.qodon.optimizers.classical_ga import GeneticAlgorithm
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-co", "GA", "-n", "4", "-c", "1000", "-sd", "2546345746583"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = GeneticAlgorithm(parser)
            opt2 = GeneticAlgorithm(parser)
            self.assertEqual(opt.final_codons, opt2.final_codons)

    def test_GA_OneIteration(self):
        from src.qodon.optimizers.classical_ga import GeneticAlgorithm
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-co", "GA", "-n", "4", "-c", "1", "-sd", "1"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = GeneticAlgorithm(parser)
            opt2 = GeneticAlgorithm(parser)
            self.assertEqual(opt.optimization_process['sequences'][-1], opt2.optimization_process['sequences'][-1])

    def test_GA_LongerSequence(self):
        from src.qodon.optimizers.classical_ga import GeneticAlgorithm
        testargs = ["design.py", "-i", "examples/spike_trim_20.fasta", "-co", "GA", "-n", "4", "-c", "1000", "-sd", "1"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = GeneticAlgorithm(parser)
            opt2 = GeneticAlgorithm(parser)
            self.assertEqual(opt.final_codons, opt2.final_codons)

    def test_RAND_DefaultSeed(self):
        from src.qodon.optimizers.random_optimizer import RandomOptimizer
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-co", "RAND", "-n", "4", "-c", "1000", "-sd", "1"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = RandomOptimizer(parser)
            opt2 = RandomOptimizer(parser)
            self.assertEqual(opt.final_codons, opt2.final_codons)

    def test_RAND_NewSeed(self):
        from src.qodon.optimizers.random_optimizer import RandomOptimizer
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-co", "RAND", "-n", "4", "-c", "1000", "-sd", "2546345746583"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = RandomOptimizer(parser)
            opt2 = RandomOptimizer(parser)
            self.assertEqual(opt.final_codons, opt2.final_codons)

    def test_RAND_OneIteration(self):
        from src.qodon.optimizers.random_optimizer import RandomOptimizer
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-co", "RAND", "-n", "4", "-c", "1", "-sd", "1"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = RandomOptimizer(parser)
            opt2 = RandomOptimizer(parser)
            self.assertEqual(opt.optimization_process['sequences'][-1], opt2.optimization_process['sequences'][-1])

    def test_RAND_LongerSequence(self):
        from src.qodon.optimizers.random_optimizer import RandomOptimizer
        testargs = ["design.py", "-i", "examples/spike_trim_20.fasta", "-co", "RAND", "-n", "4", "-c", "1000", "-sd", "1"]
        with patch.object(sys, 'argv', testargs):
            parser = Parser()
            opt = RandomOptimizer(parser)
            opt2 = RandomOptimizer(parser)
            self.assertEqual(opt.final_codons, opt2.final_codons)

if __name__ == '__main__':
    unittest.main()
