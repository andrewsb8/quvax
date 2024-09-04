import pytest
from src.config.fold_config import FoldConfig


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
