import pytest


def test_NoTranslate(mock_optimizer, caplog):
    """
    Test error raised when final codon sequence does not translate to input protein sequence

    """
    mock_optimizer.final_codons = "ATGTTCGTATTCTTAGTGTTACTGCCGCTCGTA"
    mock_optimizer._verify_dna(mock_optimizer.final_codons)
    log_entry = (
        "src.logging.logging",
        40,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Error: Codon sequence did not translate properly! Sequence: ATGTTCGTATTCTTAGTGTTACTGCCGCTCGTA",
    )
    assert log_entry in caplog.record_tuples


def test_WillTranslate(mock_optimizer, caplog):
    """
    Test to verify correct codon sequence will correctly translate to input protein sequence

    """
    mock_optimizer.final_codons = "GGCGGCGGGAAC"
    mock_optimizer._verify_dna(mock_optimizer.final_codons)
    log_entry = (
        "src.logging.logging",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Codon sequence translated properly.",
    )
    assert log_entry in caplog.record_tuples
