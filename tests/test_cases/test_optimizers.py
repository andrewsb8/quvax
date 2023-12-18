import sys
import os
import pytest
from src.params.parser import Parser

def test_TFDE():
    from src.qodon.optimizers.tf_differential_evo import TfDiffEv
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-n", "4", "-c", "4", "-ms", "2", "-co", "TFDE"]
    parser = Parser(testargs)
    parser.log.info("test_TFDE")
    opt = TfDiffEv(parser)
    with pytest.logs(level='INFO'):
        opt._verify_dna(opt.final_codons)

def test_RAND():
    from src.qodon.optimizers.random_optimizer import RandomOptimizer
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-n", "4", "-c", "4", "-ms", "2", "-co", "RAND"]
    parser = Parser(testargs)
    parser.log.info("test_RAND")
    opt = RandomOptimizer(parser)
    with pytest.logs(level='INFO'):
        opt._verify_dna(opt.final_codons)

def test_GA():
    from src.qodon.optimizers.classical_ga import GeneticAlgorithm
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-n", "4", "-c", "4", "-ms", "2", "-co", "GA"]
    parser = Parser(testargs)
    parser.log.info("test_GA")
    opt = GeneticAlgorithm(parser)
    with pytest.logs(level='INFO'):
        opt._verify_dna(opt.final_codons)
