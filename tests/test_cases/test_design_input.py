import pytest
from src.config.design_config import DesignConfig
from src.exceptions.exceptions import InvalidSequenceError


def test_input_seq_str():
    """
    Test if _validate() will throw error for a numeric input string

    """
    testargs = ["-i", "tests/test_files/test_sequences/integer_sequence.fasta"]
    with pytest.raises(InvalidSequenceError):
        DesignConfig(testargs)


def test_input_seq_partial_str():
    """
    Test if _validate() will throw error for an alphanumeric input string

    """
    testargs = ["-i", "tests/test_files/test_sequences/some_integers.fasta"]
    with pytest.raises(InvalidSequenceError):
        DesignConfig(testargs)


def test_input_seq_wrong_letter():
    """
    Test if _validate() will throw error for input string with letters not
    assigned to an amino acid

    """
    testargs = ["-i", "tests/test_files/test_sequences/non_aminoacid_letter.fasta"]
    with pytest.raises(InvalidSequenceError):
        DesignConfig(testargs)


def test_input_seq_warning(caplog):
    """
    Test if _validate() will produce a warning when it detects a sequence
    which looks like DNA instead of amino acids

    """

    testargs = ["-i", "tests/test_files/test_sequences/GAG.fasta"]
    DesignConfig(testargs)
    warning_entry = (
        "src.logging.logging",
        30,  # 30 indicates WARNING, 20 indicates INFO
        "Input protein sequence looks like an DNA sequence!",
    )
    assert warning_entry in caplog.record_tuples


def test_target_start_codon():
    """
    Test if _validate() will produce an error if start codon is in target and M is not in protein sequence

    """

    testargs = ["-i", "tests/test_files/test_sequences/GAG.fasta", "-t", "AUGAAAAAA"]
    with pytest.raises(ValueError):
        DesignConfig(testargs)


def test_target_start_codon_withM():
    """
    Test if _validate() will produce an error if there are more start codons than amino acids M

    """

    testargs = ["-i", "tests/test_files/test_sequences/GMG.fasta", "-t", "AUGAUGAAA"]
    with pytest.raises(ValueError):
        DesignConfig(testargs)


def test_target_start_codon_withM_noerror(caplog):
    """
    Test if _validate() will not produce an error if there are an equal number of start codons and amino acids M.
    Pass condition is observation of log output in method after _validate()

    """

    testargs = ["-i", "tests/test_files/test_sequences/GMG.fasta", "-t", "AAAAUGAAA"]
    DesignConfig(testargs)
    log_entry = (
        "src.logging.logging",
        10,  # 30 indicates WARNING, 20 indicates INFO, 10 DEBUG
        "\n\nList of Parameters:",
    )
    assert log_entry in caplog.record_tuples


def test_target_stop_codon():
    """
    Test if _validate() will produce an error if stop codons are in the target sequence

    """

    stop_codons = ["UAAAAAAAA", "UAGAAAAAA", "UGAAAAAAA"]
    for stop in stop_codons:
        testargs = ["-i", "tests/test_files/test_sequences/GAG.fasta", "-t", stop]
        with pytest.raises(ValueError):
            DesignConfig(testargs)
