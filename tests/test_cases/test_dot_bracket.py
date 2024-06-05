import pytest
from tests.conftest import MockFolder


def test_hairpin_dot_bracket():
    """
    Test that correct dot-bracket secondary structure works

    Test case taken from here: https://www.researchgate.net/publication/362779104_RNA_secondary_structure_factorization_in_prime_tangles/figures

    """
    sequence = "AUAGCUCAGUGGUAGAGCA"
    ref_sec_struct = "...((((.......))))."

    #manually assign stem
    stem = [(4, 18, 4)]
    folder = MockFolder(None)
    folder._stems_to_dot_bracket(19, stem)
    assert ref_sec_struct == folder.dot_bracket
