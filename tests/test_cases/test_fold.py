import pytest
from src.config.fold_config import FoldConfig


def test_scoring_func():
    """
    Test the scoring function on an input RNA structure.
    This is a long test and some refactoring should probably occur
    to shorten its length.

    """
    from src.rna_folding.rna_folder import RNAFolder

    # NOTE: this sequence won't be used to calculate energy.
    #       it just is used to create the config which contain
    #       default parameters for the hamiltonian or scoring
    #       function. Don't want to rely on random seeds to
    #       produce the same structure on each machine for
    #       the test to pass
    testargs = [
        "-i",
        "AUG",
    ]
    config = FoldConfig(testargs)
    rna_fold_obj = RNAFolder(config)

    ct_file = "tests/test_files/test_structures/trial.ct"
    structure = rna_fold_obj._ct_to_dataframe(ct_file)
    seq = rna_fold_obj._get_sequence_from_connect_table(structure)
    rna_fold_obj.stems = rna_fold_obj._connect_table_to_stems(len(seq), structure)
    rna_fold_obj.len_stem_list = len(rna_fold_obj.stems)
    idx = [i for i in range(rna_fold_obj.len_stem_list)]

    rna_fold_obj._compute_h_and_J()
    score = rna_fold_obj._calc_score(idx)

    assert score == -12.0


def test_fold_SA(caplog):
    """
    Test that fold.py will run with the Simulated Annealing Folder, including all
    outputs, with no errors.

    """
    from src.rna_folding.rna_folders.simulated_annealer import SimulatedAnnealer

    testargs = [
        "-i",
        "AUGACUAGGUAUCUAUCUUAU",
        "-r",
        "2000",
        "-ot",
        "all",
        "-o",
        "quvax",
    ]
    parser = FoldConfig(testargs)
    folder = SimulatedAnnealer(parser)
    folder._fold(folder.config.args.input, post_process=True)
    log_entry = (
        "src.logging.logging",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Folded secondary structure: ((((...((((..))))))))",
    )
    print(caplog.record_tuples)
    assert log_entry in caplog.record_tuples


def test_fold_MC(caplog):
    """
    Test that fold.py will run with the Monte Carlo Folder, including all
    outputs, with no errors.

    """
    from src.rna_folding.rna_folders.classical_mc import MC

    testargs = [
        "-i",
        "AUGACUAGGUAUCUAUCUUAU",
        "-r",
        "2000",
        "-ot",
        "all",
        "-o",
        "quvax",
    ]
    parser = FoldConfig(testargs)
    folder = MC(parser)
    folder._fold(folder.config.args.input, post_process=True)
    log_entry = (
        "src.logging.logging",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Folded secondary structure: ((((...((((..))))))))",
    )
    print(caplog.record_tuples)
    assert log_entry in caplog.record_tuples


def test_fold_ES(caplog):
    """
    Test that fold.py will run with the Exact Solver Folder, including all
    outputs, with no errors.

    """
    from src.rna_folding.rna_folders.exact_solver import ExactSolver

    testargs = [
        "-i",
        "AUGACUAGGUAUCUAUCUUAU",
        "-r",
        "2000",
        "-ot",
        "all",
        "-o",
        "quvax",
    ]
    parser = FoldConfig(testargs)
    folder = ExactSolver(parser)
    folder._fold(folder.config.args.input, post_process=True)
    log_entry = (
        "src.logging.logging",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Folded secondary structure: ((((...((((..))))))))",
    )
    print(caplog.record_tuples)
    assert log_entry in caplog.record_tuples


def test_fold_SA_no_stems(caplog):
    """
    Test that fold.py will run, including all outputs, with a sequence where
    no stems are found

    """
    from src.rna_folding.rna_folders.simulated_annealer import SimulatedAnnealer

    testargs = ["-i", "AAAAAAAAAAAAA", "-r", "2000", "-ot", "all", "-o", "quvax"]
    parser = FoldConfig(testargs)
    folder = SimulatedAnnealer(parser)
    folder._fold(folder.config.args.input, post_process=True)
    log_entry = (
        "src.logging.logging",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Folded secondary structure: .............",
    )
    print(caplog.record_tuples)
    assert log_entry in caplog.record_tuples
