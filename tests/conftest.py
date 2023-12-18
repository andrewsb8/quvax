import sys
import pytest
import unittest
from unittest.mock import patch
from src.params.parser import Parser
from src.qodon.optimizer import CodonOptimizer

class MockOptimizer(CodonOptimizer):
    def __init__(self, config):
        super().__init__(config)

    def _optimize(self):
        pass

@pytest.fixture
def mock_optimizer():
    """Mock optimizer"""
    testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-n", "5"]
    with patch.object(sys, 'argv', testargs):
        parser = Parser()
        opt = MockOptimizer(parser)
    return opt
