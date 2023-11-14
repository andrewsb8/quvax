import os
import sys
import logging
from src.params.parser import Parser

if __name__ == "__main__":
    config = Parser()
    config.log = logging.getLogger(__name__)
    config.log.setLevel(logging.DEBUG)
    handler = logging.FileHandler(config.args.log_file_name, mode='w+')
    config.log.addHandler(handler)

    if config.args.codon_optimizer == "TFDE":
        from src.qodon.optimizers.tf_differential_evo import TfDiffEv
        TfDiffEv(config)
    elif config.args.codon_optimizer == "GA":
        from src.qodon.optimizers.classical_ga import GeneticAlgorithm
        GeneticAlgorithm(config)
    elif config.args.codon_optimizer == "RAND":
        from src.qodon.optimizers.random_optimizer import RandomOptimizer
        RandomOptimizer(config)
    else:
        print("Error here.")
