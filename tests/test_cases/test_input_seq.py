import sys
import os
import unittest
from unittest.mock import patch
from src.params.parser import Parser
from Bio import SeqIO

class TestInputSeq(unittest.TestCase):
    """
    Test properties of the input sequence

    """
    def test_input_seq_str(self):
        """
        Test if input sequence is a string

        """
        testargs = ["design.py", "-i", "tests/test_sequences/integer_sequence.fasta"]
        with patch.object(sys, 'argv', testargs):
            with self.assertRaises(TypeError):
                Parser()

if __name__ == '__main__':
    unittest.main()
