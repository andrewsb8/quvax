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


def test_compare_identical_ct():
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
        analysis.metrics.sensitivity == 1
        and analysis.metrics.specificity == 1
        and analysis.metrics.f1 == 1
        and analysis.metrics.pos_predict_val == 1
    )


def test_compare_diff_ct():
    """
    Test compare connectivity tables for the different structures

    """
    from src.analysis.analyses.compare_ct import CompareCT

    testargs = [
        "compare_ct",
        "-i",
        "tests/test_files/test_analysis/test_compare_ct/test_connect_table_diff.ct",
        "-r",
        "tests/test_files/test_analysis/test_compare_ct/test_connect_table.ct",
    ]
    config = AnalysisConfig(testargs)
    analysis = CompareCT(config)

    assert (
        analysis.metrics.sensitivity == 0.9
        and round(analysis.metrics.specificity, 2) == 0.33
        and analysis.metrics.f1 == 0.9
        and analysis.metrics.pos_predict_val == 0.9
    )


def test_compare_nostem_ct():
    """
    Test compare connectivity tables where one table has
    no stems (no true or false positives)

    """
    from src.analysis.analyses.compare_ct import CompareCT

    testargs = [
        "compare_ct",
        "-i",
        "tests/test_files/test_analysis/test_compare_ct/test_connect_table_nostems.ct",
        "-r",
        "tests/test_files/test_analysis/test_compare_ct/test_connect_table.ct",
    ]
    config = AnalysisConfig(testargs)
    analysis = CompareCT(config)

    assert (
        analysis.metrics.sensitivity == 0
        and analysis.metrics.specificity == 1
        and analysis.metrics.f1 == 0
        and analysis.metrics.pos_predict_val == None
    )


def test_compare_allstem_ct():
    """
    Test compare connectivity tables where both structures are the same
    and all pairs are in stems (no true or false negatives)

    """
    from src.analysis.analyses.compare_ct import CompareCT

    testargs = [
        "compare_ct",
        "-i",
        "tests/test_files/test_analysis/test_compare_ct/test_connect_table_allstems.ct",
        "-r",
        "tests/test_files/test_analysis/test_compare_ct/test_connect_table_allstems.ct",
    ]
    config = AnalysisConfig(testargs)
    analysis = CompareCT(config)

    assert (
        analysis.metrics.sensitivity == 1
        and analysis.metrics.specificity == None
        and analysis.metrics.f1 == 1
        and analysis.metrics.pos_predict_val == 1
    )


def test_base_pair_ranges():
    """
    Test output for analysis computing average, max, min base pair
    ranges, counted in number of bases

    """
    from src.analysis.analyses.base_pair_ranges import BasePairRanges

    testargs = [
        "base_pair_ranges",
        "-i",
        "tests/test_files/test_structures/trial.ct",
    ]
    config = AnalysisConfig(testargs)
    analysis = BasePairRanges(config)

    assert (
        analysis.avg_range == 4.5
        and analysis.min_range == 2
        and analysis.max_range == 7
        and analysis.seq_len == 11
    )


def test_compute_energy():
    """
    Test output for analysis to compute energy of input structure

    """
    from src.analysis.analyses.compute_energy import ComputeEnergy

    testargs = [
        "compute_energy",
        "-i",
        "tests/test_files/test_structures/trial.ct",
        "-ms",
        "2",
        "-ml",
        "2",
    ]
    config = AnalysisConfig(testargs)
    analysis = ComputeEnergy(config)
    assert analysis.score == -12


def test_k_neighbor_energy():
    """
    Test output for analysis where energies of structure
    "neighbors" (perturbing system by adding/removing one
    stem)

    """
    from src.analysis.analyses.k_neighbor_energy import KNeighborEnergySearch

    testargs = [
        "k_neighbor_energy",
        "-i",
        "tests/test_files/test_structures/trial.ct",
        "-ms",
        "2",
        "-ml",
        "2",
    ]
    config = AnalysisConfig(testargs)
    analysis = KNeighborEnergySearch(config)
    assert analysis.energies == [-12.0, -4, 999984.0, 999984.0, -4]


def test_base_pair_types():
    """
    Test output for analysis computing average, max, min base pair
    ranges, counted in number of bases

    """
    from src.analysis.analyses.base_pair_types import BasePairTypes

    testargs = [
        "base_pair_types",
        "-i",
        "tests/test_files/test_structures/trial.ct",
    ]
    config = AnalysisConfig(testargs)
    analysis = BasePairTypes(config)

    assert (
        analysis.no_pair == 3
        and analysis.wc_base_pairs == 6
        and analysis.wobble_base_pairs == 2
        and analysis.nonwc_base_pairs == 0
        and analysis.pseudoknots == 4
        and analysis.num_bases == 11
    )


def test_classify_stems():
    """
    Test output for analysis for classifying stems.

    """
    from src.analysis.analyses.classify_stems import ClassifyStems

    testargs = [
        "classify_stems",
        "-i",
        "tests/test_files/test_structures/trial.ct",
    ]
    config = AnalysisConfig(testargs)
    analysis = ClassifyStems(config)

    assert (
        analysis.avg_stem == 2.0
        and analysis.min_stem == 2
        and analysis.max_stem == 2
        and analysis.seq_len == 11
        and analysis.num_stems == 2
        and analysis.pseudos == 1
        and analysis.overlaps == 0
    )


def test_contact_order():
    """
    Test computation of contact order for a given input structure

    """
    from src.analysis.analyses.contact_order import ContactOrder

    testargs = [
        "contact_order",
        "-i",
        "tests/test_files/test_structures/trial.ct",
    ]
    config = AnalysisConfig(testargs)
    analysis = ContactOrder(config)
    assert round(analysis.contact_order, 2) == 0.41


def test_unfold():
    """
    Test computation of contact order for a given input structure

    """
    from src.analysis.analyses.unfold import Unfold

    testargs = [
        "unfold",
        "-i",
        "tests/test_files/test_structures/trial.ct",
    ]
    config = AnalysisConfig(testargs)
    analysis = Unfold(config)
    assert len(analysis.stems) == 0 and analysis.energy == 0
