import pytest
from copy import deepcopy
from src.params.design_parser import DesignParser
from src.qodon.optimizers.replica_exchange_mc import REMCOptimizer
from src.qodon.optimizers.metro_optimizer import MetropolisOptimizer

def test_accept_changes():
    """
    Test that an appropriate sequence change will be accepted. This test applies
    to both METRO and REMC optimizers

    """
    from src.qodon.optimizers.metro_optimizer import MetropolisOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
    ]
    parser = DesignParser(testargs)
    opt = MetropolisOptimizer(parser)
    value = opt._check_changes(1, -1, 1, 0)
    assert value == True

def test_reject_changes():
    """
    Test that an appropriate sequence change will be rejected. This test applies
    to both METRO and REMC optimizers

    """
    from src.qodon.optimizers.metro_optimizer import MetropolisOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
    ]
    parser = DesignParser(testargs)
    opt = MetropolisOptimizer(parser)
    value = opt._check_changes(1, 1e10, 1, 0)
    assert value == False

def test_number_changes(mock_analysis):
    """
    Test the correct number of codons are changed in a sequence

    """
    from src.qodon.optimizers.metro_optimizer import MetropolisOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
    ]
    parser = DesignParser(testargs)
    opt = MetropolisOptimizer(parser)
    sequence = deepcopy(opt.initial_sequences[0])
    perturbed_sequence = opt._perturb_dna(opt.initial_sequences[0], 3)
    assert mock_analysis._codon_diff_list([opt._convert_ints_to_codons(perturbed_sequence), opt._convert_ints_to_codons(sequence)])[0] == 3
