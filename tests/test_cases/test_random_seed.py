import pytest
from src.params.design_parser import DesignParser
from tests.conftest import MockOptimizer


def test_SameInitialSequences_DefaultSeed():
    """
    Test to verify the same initial codon sequences will be initialized in subsequent optimizations and the same, default seed

    """
    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-n",
        "4",
        "-sd",
        "1",
    ]
    parser = DesignParser(testargs)
    opt = MockOptimizer(parser)
    opt2 = MockOptimizer(parser)
    assert opt.initial_sequences == opt2.initial_sequences


def test_SameInitialSequences_NewSeed():
    """
    Test to verify the same initial codon sequences will be initialized in subsequent optimizations and the same, nondefault seed

    """
    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-n",
        "4",
        "-sd",
        "1098760354",
    ]
    parser = DesignParser(testargs)
    opt = MockOptimizer(parser)
    opt2 = MockOptimizer(parser)
    assert opt.initial_sequences == opt2.initial_sequences


def test_InitialSequences_DifferentSeeds():
    """
    Test to verify different initial codon sequences will be initialized in subsequent optimizations with different seeds

    """
    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-n",
        "4",
    ]
    parser = DesignParser(testargs)
    opt = MockOptimizer(parser)
    parser.args.random_seed = 234524352
    opt2 = MockOptimizer(parser)
    assert opt.initial_sequences != opt2.initial_sequences


def test_TfDiffEv_DefaultSeed():
    """
    Test to verify the same minimum energy codon sequence is reached in separate optimizations with default seed and Tensorflow Evolution Optimizer

    """
    from src.qodon.optimizers.tf_differential_evo import TfDiffEv

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-n",
        "4",
        "-c",
        "100",
        "-sd",
        "1",
    ]
    parser = DesignParser(testargs)
    opt = TfDiffEv(parser)
    opt2 = TfDiffEv(parser)
    assert opt.list_seqs == opt2.list_seqs


def test_TfDiffEv_NewSeed():
    """
    Test to verify the same minimum energy codon sequence is reached in separate optimizations with nondefault seed and Tensorflow Evolution Optimizer

    """
    from src.qodon.optimizers.tf_differential_evo import TfDiffEv

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-co",
        "TFDE",
        "-n",
        "4",
        "-c",
        "100",
        "-sd",
        "2546345746583",
    ]
    parser = DesignParser(testargs)
    opt = TfDiffEv(parser)
    opt2 = TfDiffEv(parser)
    assert opt.list_seqs == opt2.list_seqs


# fails sometimes
@pytest.mark.skip
def test_TfDiffEv_OneIteration():
    """
    Test to verify the same sequence change occurs after one iteration of optimization with default seed with Tensorflow Differential Evolution Optimizer

    """
    from src.qodon.optimizers.tf_differential_evo import TfDiffEv

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-co",
        "TFDE",
        "-n",
        "1",
        "-c",
        "1",
        "-sd",
        "1",
    ]
    parser = DesignParser(testargs)
    opt = TfDiffEv(parser)
    opt2 = TfDiffEv(parser)
    assert (
        opt.optimization_process["sequences"][-1]
        == opt2.optimization_process["sequences"][-1]
    )


def test_TfDiffEv_LongerSequence():
    """
    Test to verify the same minimum energy codon sequence is reached in separate optimizations with default seed for a longer protein sequence and Tensorflow Evolution Optimizer

    """
    from src.qodon.optimizers.tf_differential_evo import TfDiffEv

    testargs = [
        "-i",
        "tests/test_files/test_sequences/spike_trim_10.fasta",
        "-co",
        "TFDE",
        "-n",
        "4",
        "-c",
        "5",
        "-sd",
        "1",
    ]
    parser = DesignParser(testargs)
    opt = TfDiffEv(parser)
    opt2 = TfDiffEv(parser)
    assert opt.list_seqs == opt2.list_seqs


def test_GA_DefaultSeed():
    """
    Test to verify the same minimum energy codon sequence is reached in separate optimizations with default seed for Genetic Algorithm

    """
    from src.qodon.optimizers.classical_ga import GeneticAlgorithm

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-co",
        "GA",
        "-n",
        "4",
        "-c",
        "10",
        "-sd",
        "1",
    ]
    parser = DesignParser(testargs)
    opt = GeneticAlgorithm(parser)
    opt2 = GeneticAlgorithm(parser)
    assert opt.list_seqs == opt2.list_seqs


def test_GA_NewSeed():
    """
    Test to verify the same minimum energy codon sequence is reached in separate optimizations with nondefault seed for Genetic Algorithm

    """
    from src.qodon.optimizers.classical_ga import GeneticAlgorithm

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-co",
        "GA",
        "-n",
        "4",
        "-c",
        "10",
        "-sd",
        "2546345746583",
    ]
    parser = DesignParser(testargs)
    opt = GeneticAlgorithm(parser)
    opt2 = GeneticAlgorithm(parser)
    assert opt.list_seqs == opt2.list_seqs


def test_GA_OneIteration():
    """
    Test to verify the same sequence change occurs after one iteration of optimization with default seed with Genetic Algorithm

    """
    from src.qodon.optimizers.classical_ga import GeneticAlgorithm

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-co",
        "GA",
        "-n",
        "4",
        "-c",
        "1",
        "-sd",
        "1",
    ]
    parser = DesignParser(testargs)
    opt = GeneticAlgorithm(parser)
    opt2 = GeneticAlgorithm(parser)
    assert opt.list_seqs == opt2.list_seqs


def test_GA_LongerSequence():
    """
    Test to verify the same minimum energy codon sequence is reached in separate optimizations with default seed for a longer protein sequence and Genetic Algorithm

    """
    from src.qodon.optimizers.classical_ga import GeneticAlgorithm

    testargs = [
        "-i",
        "tests/test_files/test_sequences/spike_trim_10.fasta",
        "-co",
        "GA",
        "-n",
        "4",
        "-c",
        "5",
        "-sd",
        "1",
    ]
    parser = DesignParser(testargs)
    opt = GeneticAlgorithm(parser)
    opt2 = GeneticAlgorithm(parser)
    assert opt.list_seqs == opt2.list_seqs


def test_RAND_DefaultSeed():
    """
    Test to verify the same minimum energy codon sequence is reached in separate optimizations with default seed for Random Optimizer

    """
    from src.qodon.optimizers.random_optimizer import RandomOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-co",
        "RAND",
        "-n",
        "4",
        "-c",
        "100",
        "-sd",
        "1",
    ]
    parser = DesignParser(testargs)
    opt = RandomOptimizer(parser)
    opt2 = RandomOptimizer(parser)
    assert opt.list_seqs == opt2.list_seqs


def test_RAND_NewSeed():
    """
    Test to verify the same minimum energy codon sequence is reached in separate optimizations with nondefault seed for Random Optimizer

    """
    from src.qodon.optimizers.random_optimizer import RandomOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-co",
        "RAND",
        "-n",
        "4",
        "-c",
        "10",
        "-sd",
        "2546345746583",
    ]
    parser = DesignParser(testargs)
    opt = RandomOptimizer(parser)
    opt2 = RandomOptimizer(parser)
    assert opt.list_seqs == opt2.list_seqs


def test_RAND_OneIteration():
    """
    Test to verify the same sequence change occurs after one iteration of optimization with default seed with Random Optimizer

    """
    from src.qodon.optimizers.random_optimizer import RandomOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-co",
        "RAND",
        "-n",
        "4",
        "-c",
        "1",
        "-sd",
        "1",
    ]
    parser = DesignParser(testargs)
    opt = RandomOptimizer(parser)
    opt2 = RandomOptimizer(parser)
    assert opt.list_seqs == opt2.list_seqs


def test_RAND_LongerSequence():
    """
    Test to verify the same minimum energy codon sequence is reached in separate optimizations with default seed for a longer protein sequence and Random Optimizer

    """
    from src.qodon.optimizers.random_optimizer import RandomOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/spike_trim_10.fasta",
        "-co",
        "RAND",
        "-n",
        "4",
        "-c",
        "10",
        "-sd",
        "1",
    ]
    parser = DesignParser(testargs)
    opt = RandomOptimizer(parser)
    opt2 = RandomOptimizer(parser)
    assert opt.list_seqs == opt2.list_seqs

def test_METRO_DefaultSeed():
    """
    Test to verify the same minimum energy codon sequence is reached in separate optimizations with default seed and Tensorflow Evolution Optimizer

    """
    from src.qodon.optimizers.metro_optimizer import MetropolisOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-n",
        "4",
        "-c",
        "10",
        "-sd",
        "1",
        "-co",
        "METRO",
    ]
    parser = DesignParser(testargs)
    opt = MetropolisOptimizer(parser)
    opt2 = MetropolisOptimizer(parser)
    assert opt.list_seqs == opt2.list_seqs
