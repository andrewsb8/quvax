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
        self.parser = Parser()
        self.parser.seq = str(SeqIO.read('tests/test_sequences/integer_sequence.fasta','fasta').seq)
        self.assertRaises(TypeError, parser._validate(self.parser))

if __name__ == '__main__':
    unittest.main()
