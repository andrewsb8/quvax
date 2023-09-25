import sys
import os
import unittest
from unittest.mock import patch
from src.params.parser import Parser

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
            with self.assertRaises(SystemExit):
                Parser()

    def test_input_seq_partial_str(self):
        """
        Test if _validate() will throw error for an alphanumeric input string

        """
        testargs = ["design.py", "-i", "tests/test_sequences/some_integers.fasta"]
        with patch.object(sys, 'argv', testargs):
            with self.assertRaises(SystemExit):
                Parser()

    def test_input_seq_wrong_letter(self):
        """
        Test if _validate() will throw error for input string with letters not
        assigned to an amino acid

        """
        testargs = ["design.py", "-i", "tests/test_sequences/non_aminoacid_letter.fasta"]
        with patch.object(sys, 'argv', testargs):
            with self.assertRaises(SystemExit):
                Parser()

if __name__ == '__main__':
    unittest.main()
