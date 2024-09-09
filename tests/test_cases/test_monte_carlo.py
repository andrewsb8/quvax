import pytest
from copy import deepcopy
from src.config.design_config import DesignConfig


def test_accept_changes():
    """
    Test that an appropriate sequence change will be accepted. This test applies
    to both METRO and REMC optimizers

    """
    from src.codon_opt.optimizers.metro_optimizer import MetropolisOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
    ]
    parser = DesignConfig(testargs)
    opt = MetropolisOptimizer(parser)
    value = opt._check_changes(1, -1, 1, 0)
    assert value == True


def test_accept_changes_diff_betas():
    """
    Test that an appropriate sequence change will be accepted when each energy
    is associated with a different temperature. This simulates checking if an
    exchange occurs in REMC.

    """
    from src.codon_opt.optimizers.metro_optimizer import MetropolisOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
    ]
    parser = DesignConfig(testargs)
    opt = MetropolisOptimizer(parser)
    value = opt._check_changes(1, -2, 2, -1)
    assert value == True


def test_reject_changes():
    """
    Test that an appropriate sequence change will be rejected. This test applies
    to both METRO and REMC optimizers

    """
    from src.codon_opt.optimizers.metro_optimizer import MetropolisOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
    ]
    parser = DesignConfig(testargs)
    opt = MetropolisOptimizer(parser)
    value = opt._check_changes(1, 1e10, 1, 0)
    assert value == False


def test_reject_changes_diff_betas():
    """
    Test that an appropriate sequence change will be rejected when each energy
    is associated with a different temperature. This simulates checking if an
    exchange occurs in REMC.

    """
    from src.codon_opt.optimizers.metro_optimizer import MetropolisOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
    ]
    parser = DesignConfig(testargs)
    opt = MetropolisOptimizer(parser)
    value = opt._check_changes(1, -1, 2, -2)
    assert value == False


def test_number_changes(mock_analysis):
    """
    Test the correct number of codons are changed in a sequence

    """
    from src.codon_opt.optimizers.metro_optimizer import MetropolisOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
    ]
    parser = DesignConfig(testargs)
    opt = MetropolisOptimizer(parser)
    sequence = deepcopy(opt.initial_sequences[0])
    perturbed_sequence = opt._perturb_dna(opt.initial_sequences[0], 3)
    assert (
        mock_analysis._codon_diff_list(
            [
                opt._convert_ints_to_codons(perturbed_sequence),
                opt._convert_ints_to_codons(sequence),
            ]
        )[0]
        == 3
    )


def test_random_sequence_generation():
    """
    Test the randomly generated sequence in the Metropolis function works correctly

    """
    from src.codon_opt.optimizers.metro_optimizer import MetropolisOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-p",
        "1",
        "-co",
        "METRO",
    ]
    parser = DesignConfig(testargs)
    opt = MetropolisOptimizer(parser)
    opt.seq_rejections = 1000  # trigger random sequence generation
    opt.energies = [0]
    seq_copy = deepcopy(opt.initial_sequences)
    opt._metropolis_iteration(opt.initial_sequences, opt.energies, [""])
    assert seq_copy[0] != opt.initial_sequences[0]  # confirm sequence is different


def test_accepted_exchange_even(mock_analysis):
    """
    Test exchanges occur correctly in REMC

    """
    from src.codon_opt.optimizers.replica_exchange_mc import REMCOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-p",
        "4",
    ]
    parser = DesignConfig(testargs)
    opt = REMCOptimizer(parser)
    opt.beta_list = [1, 2, 3, 4]  # higher beta -> lower temp!
    energies = [-25, -51, -150, -100]  # expect third and fourth to be exchanged
    # below doesn't actually have to be sequences, just validating elements of the
    # list are being rearranged
    members = ["first", "second", "third", "fourth"]
    opt._attempt_exchanges(2, members, energies)
    assert members == ["first", "second", "fourth", "third"]


def test_reject_exchange_even(mock_analysis):
    """
    Test exchanges occur correctly in REMC

    """
    from src.codon_opt.optimizers.replica_exchange_mc import REMCOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-p",
        "4",
    ]
    parser = DesignConfig(testargs)
    opt = REMCOptimizer(parser)
    opt.beta_list = [1, 2, 3, 4]  # higher beta -> lower temp!
    energies = [-25, -51, -100, -100]  # expect third and fourth to be exchanged
    # below doesn't actually have to be sequences, just validating elements of the
    # list are being rearranged
    members = ["first", "second", "third", "fourth"]
    opt._attempt_exchanges(2, members, energies)
    assert members == ["first", "second", "third", "fourth"]


def test_accepted_exchange_odd(mock_analysis):
    """
    Test exchanges occur correctly in REMC

    """
    from src.codon_opt.optimizers.replica_exchange_mc import REMCOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-p",
        "4",
    ]
    parser = DesignConfig(testargs)
    opt = REMCOptimizer(parser)
    opt.beta_list = [1, 2, 3, 4]  # higher beta -> lower temp!
    energies = [-25, -60, -33, -100]  # expect third and fourth to be exchanged
    # below doesn't actually have to be sequences, just validating elements of the
    # list are being rearranged
    members = ["first", "second", "third", "fourth"]
    opt._attempt_exchanges(1, members, energies)
    assert members == ["first", "third", "second", "fourth"]
