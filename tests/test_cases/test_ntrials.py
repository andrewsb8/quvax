import sys
import os
import pytest
from src.params.parser import Parser
from tests.conftest import MockOptimizer

def test_lenNTrials_5(mock_optimizer):
    assert len(mock_optimizer.initial_sequences) == 5

def test_lenNTrials_str():
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-n", "2"]
    with pytest.raises(TypeError):
        parser = Parser(testargs)
        parser.args.n_trials = str(parser.args.n_trials)
        opt = MockOptimizer(parser)
