import sys
import os
import pytest
from src.params.parser import Parser
from tests.conftest import MockOptimizer

def test_SameInitialSequences_DefaultSeed():
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-n", "4", "-sd", "1"]
    parser = Parser(testargs)
    opt = MockOptimizer(parser)
    opt2 = MockOptimizer(parser)
    assert opt.initial_sequences == opt2.initial_sequences

def test_SameInitialSequences_NewSeed():
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-n", "4", "-sd", "1098760354"]
    parser = Parser(testargs)
    opt = MockOptimizer(parser)
    opt2 = MockOptimizer(parser)
    assert opt.initial_sequences == opt2.initial_sequences

def test_TfDiffEv_DefaultSeed():
    from src.qodon.optimizers.tf_differential_evo import TfDiffEv
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-n", "4", "-c", "1000", "-sd", "1"]
    parser = Parser(testargs)
    opt = TfDiffEv(parser)
    opt2 = TfDiffEv(parser)
    assert opt.initial_sequences == opt2.initial_sequences

def test_TfDiffEv_NewSeed():
    from src.qodon.optimizers.tf_differential_evo import TfDiffEv
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-co", "TFDE", "-n", "4", "-c", "1000", "-sd", "2546345746583"]
    parser = Parser(testargs)
    opt = TfDiffEv(parser)
    opt2 = TfDiffEv(parser)
    assert opt.initial_sequences == opt2.initial_sequences

def test_TfDiffEv_OneIteration():
    from src.qodon.optimizers.tf_differential_evo import TfDiffEv
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-co", "TFDE", "-n", "4", "-c", "1", "-sd", "1"]
    parser = Parser(testargs)
    opt = TfDiffEv(parser)
    opt2 = TfDiffEv(parser)
    assert opt.optimization_process['sequences'][-1] == opt2.optimization_process['sequences'][-1]

def test_TfDiffEv_LongerSequence():
    from src.qodon.optimizers.tf_differential_evo import TfDiffEv
    testargs = ["-i", "examples/spike_trim_20.fasta", "-co", "TFDE", "-n", "4", "-c", "1000", "-sd", "1"]
    parser = Parser(testargs)
    opt = TfDiffEv(parser)
    opt2 = TfDiffEv(parser)
    assert opt.final_codons == opt2.final_codons

def test_GA_DefaultSeed():
    from src.qodon.optimizers.classical_ga import GeneticAlgorithm
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-co", "GA", "-n", "4", "-c", "1000", "-sd", "1"]
    parser = Parser(testargs)
    opt = GeneticAlgorithm(parser)
    opt2 = GeneticAlgorithm(parser)
    assert opt.final_codons == opt2.final_codons

def test_GA_NewSeed():
    from src.qodon.optimizers.classical_ga import GeneticAlgorithm
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-co", "GA", "-n", "4", "-c", "1000", "-sd", "2546345746583"]
    parser = Parser(testargs)
    opt = GeneticAlgorithm(parser)
    opt2 = GeneticAlgorithm(parser)
    assert opt.final_codons == opt2.final_codons

def test_GA_OneIteration():
    from src.qodon.optimizers.classical_ga import GeneticAlgorithm
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-co", "GA", "-n", "4", "-c", "1", "-sd", "1"]
    parser = Parser(testargs)
    opt = GeneticAlgorithm(parser)
    opt2 = GeneticAlgorithm(parser)
    assert opt.optimization_process['sequences'][-1] == opt2.optimization_process['sequences'][-1]

def test_GA_LongerSequence():
    from src.qodon.optimizers.classical_ga import GeneticAlgorithm
    testargs = ["-i", "examples/spike_trim_20.fasta", "-co", "GA", "-n", "4", "-c", "1000", "-sd", "1"]
    parser = Parser(testargs)
    opt = GeneticAlgorithm(parser)
    opt2 = GeneticAlgorithm(parser)
    assert opt.final_codons == opt2.final_codons

def test_RAND_DefaultSeed():
    from src.qodon.optimizers.random_optimizer import RandomOptimizer
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-co", "RAND", "-n", "4", "-c", "1000", "-sd", "1"]
    parser = Parser(testargs)
    opt = RandomOptimizer(parser)
    opt2 = RandomOptimizer(parser)
    assert opt.final_codons == opt2.final_codons

def test_RAND_NewSeed():
    from src.qodon.optimizers.random_optimizer import RandomOptimizer
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-co", "RAND", "-n", "4", "-c", "1000", "-sd", "2546345746583"]
    parser = Parser(testargs)
    opt = RandomOptimizer(parser)
    opt2 = RandomOptimizer(parser)
    assert opt.final_codons == opt2.final_codons

def test_RAND_OneIteration():
    from src.qodon.optimizers.random_optimizer import RandomOptimizer
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-co", "RAND", "-n", "4", "-c", "1", "-sd", "1"]
    parser = Parser(testargs)
    opt = RandomOptimizer(parser)
    opt2 = RandomOptimizer(parser)
    assert opt.optimization_process['sequences'][-1] == opt2.optimization_process['sequences'][-1]

def test_RAND_LongerSequence():
    from src.qodon.optimizers.random_optimizer import RandomOptimizer
    testargs = ["-i", "examples/spike_trim_20.fasta", "-co", "RAND", "-n", "4", "-c", "1000", "-sd", "1"]
    parser = Parser(testargs)
    opt = RandomOptimizer(parser)
    opt2 = RandomOptimizer(parser)
    assert opt.final_codons == opt2.final_codons
