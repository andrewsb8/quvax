import pytest
import sqlite3
from src.config.design_config import DesignConfig
from tests.conftest import MockOptimizer


def test_target_detected(mock_optimizer, caplog):
    """
    Test to verify the a target codon sequence will be detected in a database
    output from design.py

    """
    mock_optimizer.config.args.target = "AUGUUUGUGUUCCUAGUUUUAUUGCCCCUA"
    mock_optimizer.config.db = sqlite3.connect("tests/test_files/test_design/quvax.db")
    mock_optimizer.config.db_cursor = mock_optimizer.config.db.cursor()

    mock_optimizer._check_target()
    log_entry = (
        "src.logging.logging",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "The target codon sequence is in the list of minimum free energy sequences!",
    )
    assert log_entry in caplog.record_tuples


def test_target_detected_not_mfe(mock_optimizer, caplog):
    """
    Test to verify the a target codon sequence will be detected in a database
    output from design.py but is not in list of minimum free energy sequence

    """
    mock_optimizer.config.args.target = "AUGUUUGUGUUUCUUGUUUUGCUUCCACUU"
    mock_optimizer.config.db = sqlite3.connect("tests/test_files/test_design/quvax.db")
    mock_optimizer.config.db_cursor = mock_optimizer.config.db.cursor()

    mock_optimizer._check_target()
    log_entry = (
        "src.logging.logging",
        30,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "The target codon sequence was sampled but was not the lowest free energy sequence.",
    )
    assert log_entry in caplog.record_tuples


def test_target_not_detected(mock_optimizer, caplog):
    """
    Test to verify the a target codon sequence will be detected in a database
    output from design.py but is not in list of minimum free energy sequence

    """
    mock_optimizer.config.args.target = (
        "ATGTTCGTGTTCCTCGTTCTTCTCCCCCTTGTATCCAGTCAGTGTGTAAATCTGACCAAA"
    )
    mock_optimizer.config.db = sqlite3.connect("tests/test_files/test_design/quvax.db")
    mock_optimizer.config.db_cursor = mock_optimizer.config.db.cursor()

    mock_optimizer._check_target()
    log_entry = (
        "src.logging.logging",
        40,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "The target codon sequence was NOT sampled.",
    )
    assert log_entry in caplog.record_tuples
