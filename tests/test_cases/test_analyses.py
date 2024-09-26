import pytest
from src.config.analysis_config import AnalysisConfig
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
    config = AnalysisConfig(testargs)
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
    config = AnalysisConfig(testargs)
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
    config = AnalysisConfig(testargs)
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
    config = AnalysisConfig(testargs)
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
    config = AnalysisConfig(testargs)
    analysis = FreeEnergyGeneration(config)

    assert (
        filecmp.cmp(
            "test_fegeneration_out.txt",
            "tests/test_files/test_analysis/test_fe_generation_output/analysis_out.txt",
        )
    ) == True


def test_compare_ct():
    """
    Test compare connectivity tables for the same file

    """
    from src.analysis.analyses.compare_ct import CompareCT

    testargs = [
        "compare_ct",
        "-i",
        "tests/test_files/test_analysis/test_compare_ct/test_connect_table.ct",
        "-r",
        "tests/test_files/test_analysis/test_compare_ct/test_connect_table.ct",
    ]
    config = AnalysisConfig(testargs)
    analysis = CompareCT(config)

    assert (
        analysis.sensitivity == 1
        and analysis.specificity == 1
        and analysis.f1 == 1
        and analysis.pos_predict_val == 1
    )
