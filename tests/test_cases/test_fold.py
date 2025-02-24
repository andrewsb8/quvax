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
        "Folded secondary structure: (((...((((....)))))))",
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
        "Folded secondary structure: [[[..((((]]]....)))).",
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
        "Folded secondary structure: ((((...(((....)))))))",
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


def test_gen_stems(caplog):
    """
    Test _gen_stems in rna_folder.py that the list of output stems is valid.
    In this case, we are checking to see that loop conditions
    TODO: how to test if invalid?

    """
    from src.rna_folding.rna_folders.simulated_annealer import SimulatedAnnealer

    # get test sequence from fasta
    fasta = open("tests/test_files/test_sequences/5s_Acetobacter-sp.-1.fasta")
    seq = fasta.readlines()[-1]
    fasta.close()

    min_stem = 3
    min_loop = 3
    testargs = ["-i", seq, "-ms", str(min_stem), "-ml", str(min_loop)]
    parser = FoldConfig(testargs)
    folder = SimulatedAnnealer(parser)
    folder._declare_stem_vars(seq)
    folder._gen_stems()

    # validate stem tuples in list folder._stems
    for stem in folder.stems:
        # validate stem separation satisfies minimum loop length
        assert folder._calc_stem_separation(stem[0], stem[1], stem[2]) >= min_loop
        # then go through and validate all stems are formed from valid base pairs
        for i in range(stem[2]):
            # minus 1 below to switch from sequence location (start index 1) to location in sequence string (start index 0)
            assert (
                folder.nseq[stem[0] + i - 1],
                folder.nseq[stem[1] - i - 1],
            ) in folder.interactions
