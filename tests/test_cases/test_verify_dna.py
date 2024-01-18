import sys
import os
import pytest
from src.params.design_parser import DesignParser


def test_NoTranslate(mock_optimizer):
    """
    Test error raised when final codon sequence does not translate to input protein sequence

    """
    mock_optimizer.final_codons = "ATGTTCGTATTCTTAGTGTTACTGCCGCTCGTA"
    with pytest.raises(ValueError):
        mock_optimizer._verify_dna(mock_optimizer.final_codons)


def test_WillTranslate(mock_optimizer, caplog):
    """
    Test to verify correct codon sequence will correctly translate to input protein sequence

    """
    mock_optimizer.final_codons = "GGCGGCGGGAAC"
    mock_optimizer._verify_dna(mock_optimizer.final_codons)
    log_entry = (
        "src.params.design_parser",
        20,  # 30 indicates WARNING, 20 indicates INFO
        "Final codon sequence translated properly.",
    )
    assert log_entry in caplog.record_tuples
