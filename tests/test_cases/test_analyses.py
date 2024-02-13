import pytest
from src.params.analysis_parser import AnalysisParser
import os
import filecmp

def test_trajectory():
    """
    Test Trajectory analysis produces correct output

    """
    from src.analysis.analyses.trajectory import Trajectory

    testargs = ["-i", "tests/test_files/test_analysis/test_trajectory.qu", "-at", "trajectory", "-o", "test_traj_out.txt"]
    config = AnalysisParser(testargs)
    analysis = Trajectory(config)

    for i in range(analysis.config.data["generation_size"]):
        assert filecmp.cmp("test_traj_out.txt-" + str(i), "tests/test_files/test_analysis/test_trajectory_output/analysis_out.txt-" + str(i)) == True
        os.remove("test_traj_out.txt-" + str(i))
