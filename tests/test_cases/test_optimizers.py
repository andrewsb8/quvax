import pytest
from src.config.design_config import DesignConfig


def test_TFDE(caplog):
    """
    Test Tensorflow Differential Evolution Optimizer runs correctly

    """
    from src.codon_opt.optimizers.tf_differential_evo import TfDiffEv

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-p",
        "4",
        "-c",
        "4",
        "-ms",
        "2",
        "-co",
        "TFDE",
    ]
    parser = DesignConfig(testargs)
    parser.log.info("test_TFDE")
    TfDiffEv(parser)._optimize()
    log_entry = (
        "src.logging.logging",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Finished parsing optimized sequences.",
    )
    assert log_entry in caplog.record_tuples


def test_RAND(caplog):
    """
    Test Random Optimizer runs correctly

    """
    from src.codon_opt.optimizers.random_optimizer import RandomOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-p",
        "4",
        "-c",
        "4",
        "-ms",
        "2",
        "-co",
        "RAND",
    ]
    parser = DesignConfig(testargs)
    parser.log.info("test_RAND")
    RandomOptimizer(parser)._optimize()
    log_entry = (
        "src.logging.logging",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Finished parsing optimized sequences.",
    )
    assert log_entry in caplog.record_tuples


def test_GA(caplog):
    """
    Test Genetic Algorithm Optimizer runs correctly

    """
    from src.codon_opt.optimizers.classical_ga import GeneticAlgorithm

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-p",
        "4",
        "-c",
        "4",
        "-ms",
        "2",
        "-co",
        "GA",
    ]
    parser = DesignConfig(testargs)
    parser.log.info("test_GA")
    GeneticAlgorithm(parser)._optimize()
    log_entry = (
        "src.logging.logging",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Finished parsing optimized sequences.",
    )
    assert log_entry in caplog.record_tuples


def test_Metro(caplog):
    """
    Test Metropolis Algorithm Optimizer runs correctly

    """
    from src.codon_opt.optimizers.metro_optimizer import MetropolisOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-p",
        "4",
        "-c",
        "4",
        "-ms",
        "2",
        "-co",
        "METRO",
    ]
    parser = DesignConfig(testargs)
    parser.log.info("test_METRO")
    MetropolisOptimizer(parser)._optimize()
    log_entry = (
        "src.logging.logging",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Finished parsing optimized sequences.",
    )
    assert log_entry in caplog.record_tuples


def test_REMC(caplog):
    """
    Test Replica Exchange Monte Carlo Optimizer runs correctly

    """
    from src.codon_opt.optimizers.replica_exchange_mc import REMCOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-p",
        "4",
        "-c",
        "4",
        "-ms",
        "2",
        "-co",
        "REMC",
    ]
    parser = DesignConfig(testargs)
    parser.log.info("test_METRO")
    REMCOptimizer(parser)._optimize()
    log_entry = (
        "src.logging.logging",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Finished parsing optimized sequences.",
    )
    assert log_entry in caplog.record_tuples


def test_GA_MCSolver(caplog):
    """
    Test Genetic Algorithm Optimizer runs correctly with Monte Carlo RNA folding

    """
    from src.codon_opt.optimizers.classical_ga import GeneticAlgorithm

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-p",
        "4",
        "-c",
        "4",
        "-ms",
        "2",
        "-co",
        "GA",
        "-s",
        "MC",
    ]
    parser = DesignConfig(testargs)
    parser.log.info("test_GA")
    GeneticAlgorithm(parser)._optimize()
    log_entry = (
        "src.logging.logging",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Finished parsing optimized sequences.",
    )
    assert log_entry in caplog.record_tuples


def test_GA_ExactSolver(caplog):
    """
    Test Genetic Algorithm Optimizer runs correctly with ExactSolver for RNA
    folding

    """
    from src.codon_opt.optimizers.classical_ga import GeneticAlgorithm

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-p",
        "4",
        "-c",
        "4",
        "-ms",
        "2",
        "-co",
        "GA",
        "-s",
        "ES",
    ]
    parser = DesignConfig(testargs)
    parser.log.info("test_GA")
    GeneticAlgorithm(parser)._optimize()
    log_entry = (
        "src.logging.logging",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Finished parsing optimized sequences.",
    )
    assert log_entry in caplog.record_tuples


def test_convergence(caplog):
    """
    Test that convergence will be achieved

    """
    from src.codon_opt.optimizers.random_optimizer import RandomOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-p",
        "1",
        "-c",
        "10",
        "-co",
        "RAND",
        "-cc",
        "1",
    ]
    parser = DesignConfig(testargs)
    with pytest.raises(SystemExit) as s:
        RandomOptimizer(parser)._optimize()
    assert s.type == SystemExit
    assert s.value.code == 1
