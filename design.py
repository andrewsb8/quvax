import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' #removes tensorflow warnings
from src.params.parser import Parser
from src.qodon.optimizers.tf_differential_evo import TfDiffEv

if __name__ == "__main__":
    config = Parser()

    if config.args.codon_optimizer == "TFDE":
        TfDiffEv(config)
    else:
        print("Error here.")
