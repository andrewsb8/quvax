import pytest
from src.params.analysis_parser import AnalysisParser


def test_incorrect_input():
    """
    Test if _validate() will throw error when the input data is incorrectly formatted

    """
    testargs = [
        "-i",
        "tests/test_files/test_analysis/test_bad_trajectory.qu",
        "-at",
        "fe_landscape",
    ]
    with pytest.raises(ValueError):
        AnalysisParser(testargs)
