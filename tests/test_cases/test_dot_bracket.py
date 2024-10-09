import pytest
from tests.conftest import MockFolder


def test_no_stem_dot_bracket():
    """
    Test that the dot bracket function will produce a list of dots if no stems
    are provided

    """
    sequence = "AUAGCUCAGUGGUAGAGCA"
    ref_sec_struct = "..................."

    stem = []
    folder = MockFolder(None)
    dot_bracket = folder._stems_to_dot_bracket(len(sequence), stem)
    assert ref_sec_struct == dot_bracket


def test_hairpin_dot_bracket():
    """
    Test that correct dot-bracket secondary structure works

    Test case taken from here: https://www.researchgate.net/publication/362779104_RNA_secondary_structure_factorization_in_prime_tangles/figures

    """
    sequence = "AUAGCUCAGUGGUAGAGCA"
    ref_sec_struct = "...((((.......))))."

    # manually assign stem
    stem = [(4, 18, 4)]
    folder = MockFolder(None)
    dot_bracket = folder._stems_to_dot_bracket(len(sequence), stem)
    assert ref_sec_struct == dot_bracket


def test_pseudoknot_dot_bracket():
    """
    Test that correct dot-bracket secondary structure works with pseudoknots

    Test case taken from here: https://rnavlab.utep.edu/static/PKB_files/PKB104

    """
    sequence = "UCGUUGACAAGUACGAAAUCUUGUUA"
    ref_sec_struct = ".(((.[[[[[[.)))....]]]]]]."

    # manually assign stems
    stems = [(2, 15, 3), (6, 25, 6)]
    folder = MockFolder(None)
    dot_bracket = folder._stems_to_dot_bracket(len(sequence), stems)
    assert ref_sec_struct == dot_bracket
