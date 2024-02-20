import pytest
from src.params.design_parser import DesignParser
from tests.conftest import MockOptimizer


def test_lenNTrials_5(mock_optimizer):
    """
    Test the correct number of codon sequences are generated

    """
    assert len(mock_optimizer.initial_sequences) == 5


def test_lenNTrials_str():
    """
    Test exception when non-integer is provided for -n

    """
    testargs = ["-i", "tests/test_files/test_sequences/GGGN.fasta", "-n", "2"]
    with pytest.raises(TypeError):
        parser = DesignParser(testargs)
        parser.args.n_trials = str(parser.args.n_trials)
        MockOptimizer(parser)
