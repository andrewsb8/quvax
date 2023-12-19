import sys
import os
import pytest
from src.params.parser import Parser

def test_NoTranslate(mock_optimizer):
    mock_optimizer.final_codons = "ATGTTCGTATTCTTAGTGTTACTGCCGCTCGTA"
    with pytest.raises(ValueError):
        mock_optimizer._verify_dna(mock_optimizer.final_codons)

#need to find out how pytest works with logging
def test_WillTranslate(mock_optimizer, caplog):
    mock_optimizer.final_codons = "GGCGGCGGGAAC"
    mock_optimizer._verify_dna(mock_optimizer.final_codons)
    log_entry = (
        "src.params.parser",
        20, #30 indicates WARNING, 20 indicates INFO
        "Final codon sequence translated properly."
    )
    assert log_entry in caplog.record_tuples
