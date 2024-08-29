import pytest


def test_fold_SA(caplog):
    """
    Test that fold.py will run with the Simulated Annealing Folder, including all
    outputs, with no errors.

    """
    from src.params.fold_parser import FoldParser
    from src.rna_folding.rna_folders.simulated_annealer import SimulatedAnnealer

    testargs = ["-i", "AUGACUAGGUAUCUAUCUUAU", "-r", "2000", "-ot", "all"]
    parser = FoldParser(testargs)
    folder = SimulatedAnnealer(parser)
    folder._fold(folder.config.args.input, post_process=True)
    folder2 = SimulatedAnnealer(parser)
    folder2._fold(folder2.config.args.input, post_process=True)
    log_entry = (
        "src.params.fold_parser",
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
    from src.params.fold_parser import FoldParser
    from src.rna_folding.rna_folders.classical_mc import MC

    testargs = ["-i", "AUGACUAGGUAUCUAUCUUAU", "-r", "2000", "-ot", "all"]
    parser = FoldParser(testargs)
    folder = MC(parser)
    folder._fold(folder.config.args.input, post_process=True)
    folder2 = MC(parser)
    folder2._fold(folder2.config.args.input, post_process=True)
    log_entry = (
        "src.params.fold_parser",
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
    from src.params.fold_parser import FoldParser
    from src.rna_folding.rna_folders.exact_solver import ExactSolver

    testargs = ["-i", "AUGACUAGGUAUCUAUCUUAU", "-r", "2000", "-ot", "all"]
    parser = FoldParser(testargs)
    folder = ExactSolver(parser)
    folder._fold(folder.config.args.input, post_process=True)
    folder2 = ExactSolver(parser)
    folder2._fold(folder2.config.args.input, post_process=True)
    log_entry = (
        "src.params.fold_parser",
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
    from src.params.fold_parser import FoldParser
    from src.rna_folding.rna_folders.simulated_annealer import SimulatedAnnealer

    testargs = ["-i", "AAAAAAAAAAAAA", "-r", "2000", "-ot", "all"]
    parser = FoldParser(testargs)
    folder = SimulatedAnnealer(parser)
    folder._fold(folder.config.args.input, post_process=True)
    folder2 = SimulatedAnnealer(parser)
    folder2._fold(folder2.config.args.input, post_process=True)
    log_entry = (
        "src.params.fold_parser",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Folded secondary structure: .............",
    )
    print(caplog.record_tuples)
    assert log_entry in caplog.record_tuples