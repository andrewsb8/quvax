import shutil
import pytest
from src.config.design_config import DesignConfig


def test_resume(caplog):
    """
    Test to verify successful execution of resuming optimization

    """
    from src.codon_opt.optimizers.classical_ga import GeneticAlgorithm

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
    config = DesignConfig._resume(testargs)
    GeneticAlgorithm(config)._optimize()
    log_entry = (
        "src.logging.logging",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Finished parsing optimized sequences.",
    )
    assert log_entry in caplog.record_tuples


def test_resume_hash_fail():
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
        DesignConfig._resume(testargs)


def test_resume_hash(caplog):
    """
    Test to verify successful execution of resuming optimization by specifying
    hash in a database containing data from multiple optimizations

    """
    from src.codon_opt.optimizers.classical_ga import GeneticAlgorithm

    # copy database so test does not add information to test file
    # it will be deleted after test is completed
    shutil.copy("tests/test_files/test_design/quvax2.db", "quvax2.db")
    shutil.copy("tests/test_files/test_design/quvax.state", "quvax.state")
    testargs = [
        "-i",
        "quvax2.db",
        "--resume",
        "-e",
        "3",
        "-hv",
        "3de8de42d9",
    ]
    config = DesignConfig._resume(testargs)
    GeneticAlgorithm(config)._optimize()
    log_entry = (
        "src.logging.logging",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Finished parsing optimized sequences.",
    )
    assert log_entry in caplog.record_tuples


def test_resume_compare():
    """
    Test to verify that --resume will produce the same trajectory as an
    optimization that is done with a single execution of design.py

    """
    from src.codon_opt.optimizers.classical_ga import GeneticAlgorithm

    # first optimization
    testargs = [
        "-i",
        "tests/test_files/test_sequences/GAG.fasta",
        "-p",
        "2",
        "-c",
        "10",
        "-o",
        "first.db",
        "-co",
        "GA",
    ]
    config = DesignConfig(testargs)
    opt = GeneticAlgorithm(config)
    opt._optimize()

    # second optimization
    testargs2 = [
        "-i",
        "tests/test_files/test_sequences/GAG.fasta",
        "-p",
        "2",
        "-c",
        "5",
        "-o",
        "second.db",
        "-co",
        "GA",
    ]
    config2 = DesignConfig(testargs2)
    opt2 = GeneticAlgorithm(config2)
    opt2._optimize()

    # resume second optimization
    testargs3 = ["-i", "second.db", "--resume", "-e", "5"]
    config3 = DesignConfig._resume(testargs3)
    opt3 = GeneticAlgorithm(config3)
    opt3._optimize()

    assert opt.list_seqs == opt3.list_seqs


def test_resume_compare_MC():
    """
    Test to verify that --resume will produce the same trajectory as an
    optimization that is done with a single execution of design.py

    """
    from src.codon_opt.optimizers.classical_ga import GeneticAlgorithm

    # first optimization
    testargs = [
        "-i",
        "tests/test_files/test_sequences/GAG.fasta",
        "-p",
        "2",
        "-c",
        "10",
        "-o",
        "first.db",
        "-s",
        "MC",
        "-co",
        "GA",
    ]
    config = DesignConfig(testargs)
    opt = GeneticAlgorithm(config)
    opt._optimize()

    # second optimization
    testargs2 = [
        "-i",
        "tests/test_files/test_sequences/GAG.fasta",
        "-p",
        "2",
        "-c",
        "5",
        "-o",
        "second.db",
        "-s",
        "MC",
        "-co",
        "GA",
    ]
    config2 = DesignConfig(testargs2)
    opt2 = GeneticAlgorithm(config2)
    opt2._optimize()

    # resume second optimization
    testargs3 = ["-i", "second.db", "--resume", "-e", "5"]
    config3 = DesignConfig._resume(testargs3)
    opt3 = GeneticAlgorithm(config3)
    opt3._optimize()

    assert opt.list_seqs == opt3.list_seqs
