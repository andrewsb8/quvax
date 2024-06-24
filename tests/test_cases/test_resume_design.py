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
    shutil.copy("tests/test_files/test_design/quvax.state", "quvax.state")
    testargs = [
        "-i",
        "quvax.db",
        "--resume",
        "-e",
        "3",
    ]
    config = DesignParser._resume(testargs)
    GeneticAlgorithm(config)
    log_entry = (
        "src.params.design_parser",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Finished parsing optimized sequences.",
    )
    assert log_entry in caplog.record_tuples


def test_resume_hash_fail(caplog):
    """
    Test to verify graceful failure if specifying a hash value that does not exist

    """

    # copy database so test does not add information to test file
    # it will be deleted after test is completed
    shutil.copy("tests/test_files/test_design/quvax2.db", "quvax2.db")
    shutil.copy("tests/test_files/test_design/quvax.state", "quvax.state")
    testargs = [
        "-i",
        "quvax.db",
        "--resume",
        "-e",
        "3",
        "-hv",
        "0",
    ]
    with pytest.raises(ValueError):
        DesignParser._resume(testargs)


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
    testargs3 = ["-i", "second.db", "--resume", "-e", "5"]
    config3 = DesignParser._resume(testargs3)
    opt3 = GeneticAlgorithm(config3)

    assert opt.list_seqs == opt3.list_seqs


def test_resume_compare_MC(caplog):
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
        "-s",
        "MC",
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
        "-s",
        "MC",
    ]
    config2 = DesignParser(testargs2)
    opt2 = GeneticAlgorithm(config2)

    # resume second optimization
    testargs3 = ["-i", "second.db", "--resume", "-e", "5"]
    config3 = DesignParser._resume(testargs3)
    opt3 = GeneticAlgorithm(config3)

    assert opt.list_seqs == opt3.list_seqs
