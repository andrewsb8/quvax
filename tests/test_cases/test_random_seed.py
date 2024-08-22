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
        "-p",
        "4",
        "-sd",
        "1",
    ]
    parser = DesignParser(testargs)
    opt = MockOptimizer(parser)
    opt._optimize()
    opt2 = MockOptimizer(parser)
    opt2._optimize()
    assert opt.initial_sequences == opt2.initial_sequences


def test_SameInitialSequences_NewSeed():
    """
    Test to verify the same initial codon sequences will be initialized in subsequent optimizations and the same, nondefault seed

    """
    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-p",
        "4",
        "-sd",
        "1098760354",
    ]
    parser = DesignParser(testargs)
    opt = MockOptimizer(parser)
    opt._optimize()
    opt2 = MockOptimizer(parser)
    opt2._optimize()
    assert opt.initial_sequences == opt2.initial_sequences


def test_InitialSequences_DifferentSeeds():
    """
    Test to verify different initial codon sequences will be initialized in subsequent optimizations with different seeds

    """
    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-p",
        "4",
    ]
    parser = DesignParser(testargs)
    opt = MockOptimizer(parser)
    opt._optimize()
    parser.args.random_seed = 234524352
    opt2 = MockOptimizer(parser)
    opt2._optimize()
    assert opt.initial_sequences != opt2.initial_sequences


def test_TfDiffEv_DefaultSeed():
    """
    Test to verify the same minimum energy codon sequence is reached in separate optimizations with default seed and Tensorflow Evolution Optimizer

    """
    from src.qodon.optimizers.tf_differential_evo import TfDiffEv

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-p",
        "4",
        "-c",
        "100",
        "-sd",
        "1",
    ]
    parser = DesignParser(testargs)
    opt = TfDiffEv(parser)
    opt._optimize()
    opt2 = TfDiffEv(parser)
    opt2._optimize()
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
        "-p",
        "4",
        "-c",
        "100",
        "-sd",
        "2546345746583",
    ]
    parser = DesignParser(testargs)
    opt = TfDiffEv(parser)
    opt._optimize()
    opt2 = TfDiffEv(parser)
    opt2._optimize()
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
        "-p",
        "1",
        "-c",
        "1",
        "-sd",
        "1",
    ]
    parser = DesignParser(testargs)
    opt = TfDiffEv(parser)
    opt._optimize()
    opt2 = TfDiffEv(parser)
    opt2._optimize()
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
        "-p",
        "4",
        "-c",
        "5",
        "-sd",
        "1",
    ]
    parser = DesignParser(testargs)
    opt = TfDiffEv(parser)
    opt._optimize()
    opt2 = TfDiffEv(parser)
    opt2._optimize()
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
        "-p",
        "4",
        "-c",
        "10",
        "-sd",
        "1",
    ]
    parser = DesignParser(testargs)
    opt = GeneticAlgorithm(parser)
    opt._optimize()
    opt2 = GeneticAlgorithm(parser)
    opt2._optimize()
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
        "-p",
        "4",
        "-c",
        "10",
        "-sd",
        "2546345746583",
    ]
    parser = DesignParser(testargs)
    opt = GeneticAlgorithm(parser)
    opt._optimize()
    opt2 = GeneticAlgorithm(parser)
    opt2._optimize()
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
        "-p",
        "4",
        "-c",
        "1",
        "-sd",
        "1",
    ]
    parser = DesignParser(testargs)
    opt = GeneticAlgorithm(parser)
    opt._optimize()
    opt2 = GeneticAlgorithm(parser)
    opt2._optimize()
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
        "-p",
        "4",
        "-c",
        "5",
        "-sd",
        "1",
    ]
    parser = DesignParser(testargs)
    opt = GeneticAlgorithm(parser)
    opt._optimize()
    opt2 = GeneticAlgorithm(parser)
    opt2._optimize()
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
        "-p",
        "4",
        "-c",
        "100",
        "-sd",
        "1",
    ]
    parser = DesignParser(testargs)
    opt = RandomOptimizer(parser)
    opt._optimize()
    opt2 = RandomOptimizer(parser)
    opt2._optimize()
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
        "-p",
        "4",
        "-c",
        "10",
        "-sd",
        "2546345746583",
    ]
    parser = DesignParser(testargs)
    opt = RandomOptimizer(parser)
    opt._optimize()
    opt2 = RandomOptimizer(parser)
    opt2._optimize()
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
        "-p",
        "4",
        "-c",
        "1",
        "-sd",
        "1",
    ]
    parser = DesignParser(testargs)
    opt = RandomOptimizer(parser)
    opt._optimize()
    opt2 = RandomOptimizer(parser)
    opt2._optimize()
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
        "-p",
        "4",
        "-c",
        "10",
        "-sd",
        "1",
    ]
    parser = DesignParser(testargs)
    opt = RandomOptimizer(parser)
    opt._optimize()
    opt2 = RandomOptimizer(parser)
    opt2._optimize()
    assert opt.list_seqs == opt2.list_seqs


def test_METRO_DefaultSeed():
    """
    Test to verify the same minimum energy codon sequence is reached in separate optimizations with default seed and Tensorflow Evolution Optimizer

    """
    from src.qodon.optimizers.metro_optimizer import MetropolisOptimizer

    testargs = [
        "-i",
        "tests/test_files/test_sequences/GGGN.fasta",
        "-p",
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
    opt._optimize()
    opt2 = MetropolisOptimizer(parser)
    opt2._optimize()
    assert opt.list_seqs == opt2.list_seqs

def test_Fold_SA():
    """
    Test to verify the same minimum energy and secondary structure are reached
    on repeated runs of fold.py

    """
    from src.params.fold_parser import FoldParser
    from src.rna_folding.rna_folders.simulated_annealer import SimulatedAnnealer

    testargs = [
        "-i",
        "AUGACUAGGUAUCUAUCUUAU",
    ]
    parser = FoldParser(testargs)
    folder = SimulatedAnnealer(parser)
    folder._fold(folder.config.args.input)
    folder2 = SimulatedAnnealer(parser)
    folder2._fold(folder2.config.args.input)
    assert folder.best_score == folder2.best_score and folder.dot_bracket == folder2.dot_bracket
