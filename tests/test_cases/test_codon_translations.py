import pytest


def test_codon_translation(mock_optimizer):
    """
    Test that conversion between indices and codon strings is correct

    """
    sequence = mock_optimizer.initial_sequences[0]
    codon_string = mock_optimizer._convert_ints_to_codons(sequence)
    sequence2 = mock_optimizer._convert_codons_to_ints(codon_string)
    assert sequence == sequence2
