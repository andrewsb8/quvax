import pytest
from src.params.design_parser import DesignParser


def test_arg_codon_iterations():
    """
    Test if _validate() will throw error when -c is less than 1

    """
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-c", "0"]
    with pytest.raises(ValueError):
        DesignParser(testargs)


def test_arg_rna_iterations():
    """
    Test if _validate() will throw error when -r is less than 1

    """
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-r", "0"]
    with pytest.raises(ValueError):
        DesignParser(testargs)


def test_arg_ntrials():
    """
    Test if _validate() will throw error when -n is less than 1

    """
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-n", "0"]
    with pytest.raises(ValueError):
        DesignParser(testargs)
