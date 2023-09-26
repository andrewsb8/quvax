import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' #removes tensorflow warnings
from src.params.parser import Parser
from src.qodon.optimizers.tf_differential_evo import TfDiffEv
from src.qodon.optimizers.classical_ga import GeneticAlgorithm

if __name__ == "__main__":
    config = Parser()

    if config.args.codon_optimizer == "TFDE":
        TfDiffEv(config)
    elif config.args.codon_optimizer == "GA":
        GeneticAlgorithm(config)
    else:
        print("Error here.")
