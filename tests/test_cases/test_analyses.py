import pytest
from src.params.analysis_parser import AnalysisParser
import os
import filecmp


def test_fe_trajectory():
    """
    Test FETrajectory analysis produces correct output

    """
    from src.analysis.analyses.fe_trajectory import FreeEnergyTrajectory

    testargs = [
        "fe_trajectory",
        "-i",
        "tests/test_files/test_analysis/quvax.db",
        "-o",
        "test_traj_out.txt",
    ]
    config = AnalysisParser(testargs)
    analysis = FreeEnergyTrajectory(config)

    for i in range(analysis.config.sim_details["generation_size"]):
        assert (
            filecmp.cmp(
                "test_traj_out.txt-" + str(i),
                "tests/test_files/test_analysis/test_fe_trajectory_output/analysis_out.txt-"
                + str(i),
            )
            == True
        )
        os.remove("test_traj_out.txt-" + str(i))


def test_fe_landscape():
    """
    Test Free Energy Landscape, with respect to min energy sequence, analysis
    produces correct output

    """
    from src.analysis.analyses.fe_landscape import FreeEnergyLandscape

    testargs = [
        "fe_landscape",
        "-i",
        "tests/test_files/test_analysis/quvax.db",
        "-o",
        "test_felandscape_out.txt",
    ]
    config = AnalysisParser(testargs)
    analysis = FreeEnergyLandscape(config)

    assert (
        filecmp.cmp(
            "test_felandscape_out.txt",
            "tests/test_files/test_analysis/test_fe_landscape_output/analysis_out.txt",
        )
    ) == True


def test_codon_trajectory():
    """
    Test CodonTrajectory analysis produces correct output

    """
    from src.analysis.analyses.codon_trajectory import CodonTrajectory

    testargs = [
        "codon_trajectory",
        "-i",
        "tests/test_files/test_analysis/quvax.db",
        "-o",
        "test_traj_out.txt",
    ]
    config = AnalysisParser(testargs)
    analysis = CodonTrajectory(config)

    for i in range(analysis.config.sim_details["generation_size"]):
        assert (
            filecmp.cmp(
                "test_traj_out.txt-" + str(i),
                "tests/test_files/test_analysis/test_codon_trajectory_output/analysis_out.txt-"
                + str(i),
            )
            == True
        )
        os.remove("test_traj_out.txt-" + str(i))


def test_fe_generation():
    """
    Test Free Energy of sequences of subsequent generations analysis produces
    correct output

    """
    from src.analysis.analyses.fe_generation import FreeEnergyGeneration

    testargs = [
        "fe_generation",
        "-i",
        "tests/test_files/test_analysis/quvax.db",
        "-o",
        "test_fegeneration_out.txt",
    ]
    config = AnalysisParser(testargs)
    analysis = FreeEnergyGeneration(config)

    assert (
        filecmp.cmp(
            "test_fegeneration_out.txt",
            "tests/test_files/test_analysis/test_fe_generation_output/analysis_out.txt",
        )
    ) == True


def test_hash_analysis():
    """
    Test that the input hash value will point at the correct optimization for analysis

    """
    from src.analysis.analyses.fe_generation import FreeEnergyGeneration

    testargs = [
        "fe_generation",
        "-i",
        "tests/test_files/test_analysis/quvax2.db",
        "-o",
        "test_fegeneration_out.txt",
        "-hv",
        "7224012419568161126",
    ]
    config = AnalysisParser(testargs)
    analysis = FreeEnergyGeneration(config)

    assert (
        filecmp.cmp(
            "test_fegeneration_out.txt",
            "tests/test_files/test_analysis/test_fe_generation_output/analysis_out.txt",
        )
    ) == True
