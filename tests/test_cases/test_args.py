import sys
import os
import unittest
from unittest.mock import patch
from src.params.parser import Parser

class TestInputSeq(unittest.TestCase):
    """
    Test properties of the input sequence

    """
    def test_arg_codon_iterations(self):
        """
        Test if _validate() will throw error when -c is less than 1

        """
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-c", "0"]
        with patch.object(sys, 'argv', testargs):
            with self.assertRaises(ValueError):
                Parser()

    def test_arg_rna_iterations(self):
        """
        Test if _validate() will throw error when -r is less than 1

        """
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-r", "0"]
        with patch.object(sys, 'argv', testargs):
            with self.assertRaises(ValueError):
                Parser()

    def test_arg_ntrials(self):
        """
        Test if _validate() will throw error when -n is less than 1

        """
        testargs = ["design.py", "-i", "tests/test_sequences/GGGN.fasta", "-n", "0"]
        with patch.object(sys, 'argv', testargs):
            with self.assertRaises(ValueError):
                Parser()



if __name__ == '__main__':
    unittest.main()
