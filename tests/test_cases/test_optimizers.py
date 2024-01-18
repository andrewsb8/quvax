import sys
import os
import pytest
from src.params.design_parser import DesignParser

def test_TFDE(caplog):
    """
    Test Tensorflow Differential Evolution Optimizer runs correctly

    """
    from src.qodon.optimizers.tf_differential_evo import TfDiffEv
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-n", "4", "-c", "4", "-ms", "2", "-co", "TFDE"]
    parser = DesignParser(testargs)
    parser.log.info("test_TFDE")
    opt = TfDiffEv(parser)
    log_entry = (
        "src.params.design_parser",
        20, #30 indicates WARNING, 20 indicates INFO
        "Final codon sequence translated properly."
    )
    assert log_entry in caplog.record_tuples

def test_RAND(caplog):
    """
    Test Random Optimizer runs correctly

    """
    from src.qodon.optimizers.random_optimizer import RandomOptimizer
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-n", "4", "-c", "4", "-ms", "2", "-co", "RAND"]
    parser = DesignParser(testargs)
    parser.log.info("test_RAND")
    opt = RandomOptimizer(parser)
    log_entry = (
        "src.params.design_parser",
        20, #30 indicates WARNING, 20 indicates INFO
        "Final codon sequence translated properly."
    )
    assert log_entry in caplog.record_tuples

def test_GA(caplog):
    """
    Test Genetic Algorithm Optimizer runs correctly

    """
    from src.qodon.optimizers.classical_ga import GeneticAlgorithm
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-n", "4", "-c", "4", "-ms", "2", "-co", "GA"]
    parser = DesignParser(testargs)
    parser.log.info("test_GA")
    opt = GeneticAlgorithm(parser)
    log_entry = (
        "src.params.design_parser",
        20, #30 indicates WARNING, 20 indicates INFO
        "Final codon sequence translated properly."
    )
    assert log_entry in caplog.record_tuples
