import os
import pytest
from src.qodon.optimizer import CodonOptimizer
from src.config.design_config import DesignConfig
from src.config.fold_config import FoldConfig
from src.rna_folding.rna_folder import RNAFolder
from src.config.analysis_config import AnalysisConfig
from src.analysis.analysis import Analysis


class MockOptimizer(CodonOptimizer):
    def __init__(self, config):
        super().__init__(config)

    def _optimize(self):
        pass


@pytest.fixture
def mock_optimizer():
    """Mock optimizer"""
    testargs = ["-i", "tests/test_files/test_sequences/GGGN.fasta", "-p", "5"]
    parser = DesignConfig(testargs)
    opt = MockOptimizer(parser)
    return opt


class MockAnalysis(Analysis):
    def __init__(self, config):
        super().__init__(config)

    def _analyze(self):
        pass


@pytest.fixture
def mock_analysis():
    """
    Mock analysis object used only to access certain functions. No config is
    needed so None is passed to it.

    """
    return MockAnalysis(None)


class MockFolder(RNAFolder):
    def __init__(self, config):
        super().__init__(config)

    def _fold(self):
        pass


@pytest.fixture(autouse=True)
def test_cleanup():
    """
    Deletes files produced by tests. Exception for Trajectory which handles the
    deletion in test test_trajectory in test_analyses.py

    """
    yield
    files = [
        "quvax.db",
        "quvax2.db",
        "quvax.state",
        "first.db",
        "second.db",
        "quvax.log",
        "test_felandscape_out.txt",
        "test_fegeneration_out.txt",
        "quvax.dot",
        "quvax.ct",
        "test_traj_out.txt",
    ]
    for file in files:
        if os.path.exists(file):
            os.remove(file)
