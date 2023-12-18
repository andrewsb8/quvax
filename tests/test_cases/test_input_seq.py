import sys
import os
import logging
import pytest
from src.params.parser import Parser
from src.exceptions.exceptions import InvalidSequenceError
import warnings

def test_input_seq_str():
    """
    Test if _validate() will throw error for a numeric input string

    """
    testargs = ["-i", "tests/test_sequences/integer_sequence.fasta"]
    with pytest.raises(InvalidSequenceError):
        Parser(testargs)

def test_input_seq_partial_str():
    """
    Test if _validate() will throw error for an alphanumeric input string

    """
    testargs = ["-i", "tests/test_sequences/some_integers.fasta"]
    with pytest.raises(InvalidSequenceError):
        Parser(testargs)

def test_input_seq_wrong_letter():
    """
    Test if _validate() will throw error for input string with letters not
    assigned to an amino acid

    """
    testargs = ["-i", "tests/test_sequences/non_aminoacid_letter.fasta"]
    with pytest.raises(InvalidSequenceError):
        Parser(testargs)

#won't pass due to reference to "self". need to find out how pytest works with logging
def test_input_seq_warning():
    """
    Test if _validate() will produce a warning when it detects a sequence
    which looks like DNA instead of amino acids

    """
    testargs = ["-i", "tests/test_sequences/GAG.fasta"]
    with pytest.logs(level='WARNING'):
        Parser(testargs)
