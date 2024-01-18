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
    testargs = ["-i", "tests/test_sequences/GGGN.fasta", "-n", "5"]
    parser = DesignParser(testargs)
    opt = MockOptimizer(parser)
    return opt
