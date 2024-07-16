import pytest
from src.params.design_parser import DesignParser
from tests.conftest import MockOptimizer


def test_arg_codon_iterations():
    """
    Test if _validate() will throw error when -c is less than 1

    """
    testargs = ["-i", "tests/test_files/test_sequences/GGGN.fasta", "-c", "0"]
    with pytest.raises(ValueError):
        DesignParser(testargs)


def test_arg_rna_iterations():
    """
    Test if _validate() will throw error when -r is less than 1

    """
    testargs = ["-i", "tests/test_files/test_sequences/GGGN.fasta", "-r", "0"]
    with pytest.raises(ValueError):
        DesignParser(testargs)


def test_arg_ntrials():
    """
    Test if _validate() will throw error when -p is less than 1

    """
    testargs = ["-i", "tests/test_files/test_sequences/GGGN.fasta", "-p", "0"]
    with pytest.raises(ValueError):
        DesignParser(testargs)


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
        parser = DesignParser(testargs)
        parser.args.population_size = str(parser.args.population_size)
        MockOptimizer(parser)
