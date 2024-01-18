import pytest
from src.params.design_parser import DesignParser
from src.exceptions.exceptions import InvalidSequenceError


def test_input_seq_str():
    """
    Test if _validate() will throw error for a numeric input string

    """
    testargs = ["-i", "tests/test_sequences/integer_sequence.fasta"]
    with pytest.raises(InvalidSequenceError):
        DesignParser(testargs)


def test_input_seq_partial_str():
    """
    Test if _validate() will throw error for an alphanumeric input string

    """
    testargs = ["-i", "tests/test_sequences/some_integers.fasta"]
    with pytest.raises(InvalidSequenceError):
        DesignParser(testargs)


def test_input_seq_wrong_letter():
    """
    Test if _validate() will throw error for input string with letters not
    assigned to an amino acid

    """
    testargs = ["-i", "tests/test_sequences/non_aminoacid_letter.fasta"]
    with pytest.raises(InvalidSequenceError):
        DesignParser(testargs)


def test_input_seq_warning(caplog):
    """
    Test if _validate() will produce a warning when it detects a sequence
    which looks like DNA instead of amino acids

    """

    testargs = ["-i", "tests/test_sequences/GAG.fasta"]
    DesignParser(testargs)
    warning_entry = (
        "src.params.design_parser",
        30,  # 30 indicates WARNING, 20 indicates INFO
        "Input protein sequence looks like an DNA sequence!",
    )
    assert warning_entry in caplog.record_tuples
