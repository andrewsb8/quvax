import sys
import os
import pytest
from src.params.parser import Parser

def test_NoTranslate(mock_optimizer):
    mock_optimizer.final_codons = "ATGTTCGTATTCTTAGTGTTACTGCCGCTCGTA"
    with pytest.raises(ValueError):
        mock_optimizer._verify_dna(mock_optimizer.final_codons)

#need to find out how pytest works with logging
def test_WillTranslate(mock_optimizer):
    mock_optimizer.final_codons = "GGCGGCGGGAAC"
    with pytest.logs(level='INFO'):
        mock_optimizer._verify_dna(mock_optimizer.final_codons)
