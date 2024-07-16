import pytest
from src.params.design_parser import DesignParser


def test_TFDE(caplog):
    """
    Test Tensorflow Differential Evolution Optimizer runs correctly

    """
    from src.qodon.optimizers.tf_differential_evo import TfDiffEv

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
    parser = DesignParser(testargs)
    parser.log.info("test_TFDE")
    TfDiffEv(parser)
    log_entry = (
        "src.params.design_parser",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Finished parsing optimized sequences.",
    )
    assert log_entry in caplog.record_tuples


def test_RAND(caplog):
    """
    Test Random Optimizer runs correctly

    """
    from src.qodon.optimizers.random_optimizer import RandomOptimizer

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
    parser = DesignParser(testargs)
    parser.log.info("test_RAND")
    RandomOptimizer(parser)
    log_entry = (
        "src.params.design_parser",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Finished parsing optimized sequences.",
    )
    assert log_entry in caplog.record_tuples


def test_GA(caplog):
    """
    Test Genetic Algorithm Optimizer runs correctly

    """
    from src.qodon.optimizers.classical_ga import GeneticAlgorithm

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
    parser = DesignParser(testargs)
    parser.log.info("test_GA")
    GeneticAlgorithm(parser)
    log_entry = (
        "src.params.design_parser",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Finished parsing optimized sequences.",
    )
    assert log_entry in caplog.record_tuples


def test_Metro(caplog):
    """
    Test Metropolis Algorithm Optimizer runs correctly

    """
    from src.qodon.optimizers.metro_optimizer import MetropolisOptimizer

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
    parser = DesignParser(testargs)
    parser.log.info("test_METRO")
    MetropolisOptimizer(parser)
    log_entry = (
        "src.params.design_parser",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Finished parsing optimized sequences.",
    )
    assert log_entry in caplog.record_tuples


def test_GA_MCSolver(caplog):
    """
    Test Genetic Algorithm Optimizer runs correctly with Monte Carlo RNA folding

    """
    from src.qodon.optimizers.classical_ga import GeneticAlgorithm

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
    parser = DesignParser(testargs)
    parser.log.info("test_GA")
    GeneticAlgorithm(parser)
    log_entry = (
        "src.params.design_parser",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Finished parsing optimized sequences.",
    )
    assert log_entry in caplog.record_tuples


def test_GA_ExactSolver(caplog):
    """
    Test Genetic Algorithm Optimizer runs correctly with ExactSolver for RNA
    folding

    """
    from src.qodon.optimizers.classical_ga import GeneticAlgorithm

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
    parser = DesignParser(testargs)
    parser.log.info("test_GA")
    GeneticAlgorithm(parser)
    log_entry = (
        "src.params.design_parser",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Finished parsing optimized sequences.",
    )
    assert log_entry in caplog.record_tuples


def test_convergence(caplog):
    """
    Test Genetic Algorithm Optimizer runs correctly

    """
    from src.qodon.optimizers.random_optimizer import RandomOptimizer

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
    parser = DesignParser(testargs)
    with pytest.raises(SystemExit) as s:
        RandomOptimizer(parser)
    assert s.type == SystemExit
    assert s.value.code == 1
