import os
import pytest
from src.qodon.optimizer import CodonOptimizer
from src.params.design_parser import DesignParser


class MockOptimizer(CodonOptimizer):
    def __init__(self, config):
        super().__init__(config)

    def _optimize(self):
        pass


@pytest.fixture
def mock_optimizer():
    """Mock optimizer"""
    testargs = ["-i", "tests/test_files/test_sequences/GGGN.fasta", "-n", "5"]
    parser = DesignParser(testargs)
    opt = MockOptimizer(parser)
    return opt


@pytest.fixture(scope="session", autouse=True)
def test_cleanup():
    """
    Deletes files produced by tests. Exception for Trajectory which handles the
    deletion in test test_trajectory in test_analyses.py

    """
    yield
    files = [
        "quvax.qu",
        "quvax.log",
        "sequences.txt",
        "test_felandscape_out.txt",
        "test_fegeneration_out.txt",
    ]
    for file in files:
        if os.path.exists(file):
            os.remove(file)
