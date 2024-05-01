import shutil
import pytest
from src.params.design_parser import DesignParser


def test_resume(caplog):
    """
    Test to verify successful execution of resuming optimization

    """
    from src.qodon.optimizers.classical_ga import GeneticAlgorithm

    # copy database so test does not add information to test file
    # it will be deleted after test is completed
    shutil.copy("tests/test_files/test_design/quvax.db", "quvax.db")
    testargs = [
        "-i",
        "quvax.db",
        "--resume",
    ]
    config = DesignParser._resume(testargs)
    GeneticAlgorithm(config)
    log_entry = (
        "src.params.design_parser",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Finished parsing optimized sequences.",
    )
    assert log_entry in caplog.record_tuples


def test_resume_compare(caplog):
    """
    Test to verify that --resume will produce the same trajectory as an
    optimization that is done with a single execution of design.py

    """
    from src.qodon.optimizers.classical_ga import GeneticAlgorithm

    # first optimization
    testargs = [
        "-i",
        "tests/test_files/test_sequences/GAG.fasta",
        "-n",
        "2",
        "-c",
        "10",
        "-o",
        "first.db",
    ]
    config = DesignParser(testargs)
    opt = GeneticAlgorithm(config)

    # second optimization
    testargs2 = [
        "-i",
        "tests/test_files/test_sequences/GAG.fasta",
        "-n",
        "2",
        "-c",
        "5",
        "-o",
        "second.db",
    ]
    config2 = DesignParser(testargs2)
    opt2 = GeneticAlgorithm(config2)

    # resume second optimization
    testargs3 = ["-i", "second.db", "--resume"]
    config3 = DesignParser._resume(testargs3)
    opt3 = GeneticAlgorithm(config3)

    assert opt.n_seqs == opt3.n_seqs
