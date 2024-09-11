import pytest
from src.config.design_config import DesignConfig
from tests.conftest import MockOptimizer


def test_arg_codon_iterations():
    """
    Test if _validate() will throw error when -c is less than 1

    """
    testargs = ["-i", "tests/test_files/test_sequences/GGGN.fasta", "-c", "0"]
    with pytest.raises(ValueError):
        DesignConfig(testargs)


def test_arg_rna_iterations():
    """
    Test if _validate() will throw error when -r is less than 1

    """
    testargs = ["-i", "tests/test_files/test_sequences/GGGN.fasta", "-r", "0"]
    with pytest.raises(ValueError):
        DesignConfig(testargs)


def test_arg_ntrials():
    """
    Test if _validate() will throw error when -p is less than 1

    """
    testargs = ["-i", "tests/test_files/test_sequences/GGGN.fasta", "-p", "0"]
    with pytest.raises(ValueError):
        DesignConfig(testargs)


def test_lenNTrials_5(mock_optimizer):
    """
    Test the correct number of codon sequences are generated

    """
    assert len(mock_optimizer.initial_sequences) == 5


def test_lenNTrials_str():
    """
    Test exception when non-integer is provided for -p

    """
    testargs = ["-i", "tests/test_files/test_sequences/GGGN.fasta", "-p", "4"]
    with pytest.raises(TypeError):
        parser = DesignConfig(testargs)
        parser.args.population_size = str(parser.args.population_size)
        MockOptimizer(parser)
