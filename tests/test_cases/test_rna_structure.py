from src.rna_structure.structure_convert import StructureConvert
from src.rna_structure.structure_io import StructureIO


def test_connect_table_to_dot_bracket():
    """
    Test correct translation of connectivity tables to dot bracket notation

    """

    ct_input = "tests/test_files/test_structures/trial.ct"
    ref_dot_bracket = "(([[)).]].."

    connect_table = StructureIO()._ct_to_dataframe(ct_input)
    num_bases = connect_table["Index"].iloc[-1]
    stems = StructureConvert()._connect_table_to_stems(num_bases, connect_table)
    dot_bracket = StructureConvert()._stems_to_dot_bracket(num_bases, stems)

    assert ref_dot_bracket == dot_bracket


def test_connect_table_to_dot_bracket_long():
    """
    Test correct translation of connectivity tables to dot bracket notation
    for a sequence with a stem late in the sequence

    """

    ct_input = "tests/test_files/test_structures/trial_long.ct"
    ref_dot_bracket = "..(((((..[[[[[)))))......]]]]]"

    connect_table = StructureIO()._ct_to_dataframe(ct_input)
    num_bases = connect_table["Index"].iloc[-1]
    stems = StructureConvert()._connect_table_to_stems(num_bases, connect_table)
    dot_bracket = StructureConvert()._stems_to_dot_bracket(num_bases, stems)

    assert ref_dot_bracket == dot_bracket
