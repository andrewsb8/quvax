import os
from src.params.parser import Parser
from src.qodon.optimizers.tf_differential_evo import TfDiffEv
from src.qodon.optimizers.classical_ga import GeneticAlgorithm
from src.qodon.optimizers.random_optimizer import RandomOptimizer

if __name__ == "__main__":
    config = Parser()

    if config.args.codon_optimizer == "TFDE":
        TfDiffEv(config)
    elif config.args.codon_optimizer == "GA":
        GeneticAlgorithm(config)
    elif config.args.codon_optimizer == "RAND":
        RandomOptimizer(config)
    else:
        print("Error here.")
