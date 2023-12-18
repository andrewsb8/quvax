import sys
import os
import pytest
import unittest
from unittest.mock import patch
from src.params.parser import Parser
from tests.conftest import MockOptimizer

def test_lenNTrials_5(mock_optimizer):
    assert len(mock_optimizer.initial_sequences) == 5

@pytest.mark.parametrize()


class TestNTrials(unittest.TestCase):
    """
    Test variable n_trials will generate appropriate number of sequences

    """

    def test_lenNTrials_str(self):
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-n", "2"]
        with patch.object(sys, 'argv', testargs):
            with self.assertRaises(TypeError):
                parser = Parser()
                parser.args.n_trials = str(parser.args.n_trials)
                opt = MockOptimizer(parser)

if __name__ == '__main__':
    unittest.main()
