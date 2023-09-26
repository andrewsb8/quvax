import sys
import os
import unittest
from unittest.mock import patch
from src.params.parser import Parser
from src.params.parser import InvalidSequenceError
import warnings

class TestInputSeq(unittest.TestCase):
    """
    Test properties of the input sequence

    """
    def test_input_seq_str(self):
        """
        Test if _validate() will throw error for a numeric input string

        """
        testargs = ["design.py", "-i", "tests/test_sequences/integer_sequence.fasta"]
        with patch.object(sys, 'argv', testargs):
            with self.assertRaises(InvalidSequenceError):
                Parser()

    def test_input_seq_partial_str(self):
        """
        Test if _validate() will throw error for an alphanumeric input string

        """
        testargs = ["design.py", "-i", "tests/test_sequences/some_integers.fasta"]
        with patch.object(sys, 'argv', testargs):
            with self.assertRaises(InvalidSequenceError):
                Parser()

    def test_input_seq_wrong_letter(self):
        """
        Test if _validate() will throw error for input string with letters not
        assigned to an amino acid

        """
        testargs = ["design.py", "-i", "tests/test_sequences/non_aminoacid_letter.fasta"]
        with patch.object(sys, 'argv', testargs):
            with self.assertRaises(InvalidSequenceError):
                Parser()

    def test_input_seq_warning(self):
        """
        Test if _validate() will produce a warning when it detects a sequence
        which looks like DNA instead of amino acids

        """
        testargs = ["design.py", "-i", "tests/test_sequences/GAG.fasta"]
        with patch.object(sys, 'argv', testargs):
            with self.assertWarns(Warning):
                Parser()

if __name__ == '__main__':
    unittest.main()
