import sys
import os
import unittest
from unittest.mock import patch
from src.params.parser import Parser

class TestInputSeq(unittest.TestCase):
    """
    Test properties of the input sequence

    """
    def test_codon_iterations(self):
        """
        Test if _validate() will throw error for a numeric input string

        """
        testargs = ["design.py", "-i", "tests/test_sequences/GGG.fasta", "-c", "0"]
        with patch.object(sys, 'argv', testargs):
            with self.assertRaises(ValueError):
                Parser()

    def test_rna_iterations(self):
        """
        Test if _validate() will throw error for a numeric input string

        """
        testargs = ["design.py", "-i", "tests/test_sequences/GGG.fasta", "-r", "0"]
        with patch.object(sys, 'argv', testargs):
            with self.assertRaises(ValueError):
                Parser()

    def test_ntrials(self):
        """
        Test if _validate() will throw error for a numeric input string

        """
        testargs = ["design.py", "-i", "tests/test_sequences/GGG.fasta", "-n", "0"]
        with patch.object(sys, 'argv', testargs):
            with self.assertRaises(ValueError):
                Parser()



if __name__ == '__main__':
    unittest.main()
