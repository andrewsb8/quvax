import sys
from src.config.design_config import DesignConfig


if __name__ == "__main__":
    if "--resume" in sys.argv:
        config = DesignConfig._resume()
    else:
        config = DesignConfig()

    if config.args.codon_optimizer == "TFDE":
        from src.codon_opt.optimizers.tf_differential_evo import TfDiffEv

        TfDiffEv(config)._optimize()
    elif config.args.codon_optimizer == "GA":
        from src.codon_opt.optimizers.classical_ga import GeneticAlgorithm

        GeneticAlgorithm(config)._optimize()
    elif config.args.codon_optimizer == "RAND":
        from src.codon_opt.optimizers.random_optimizer import RandomOptimizer

        RandomOptimizer(config)._optimize()
    elif config.args.codon_optimizer == "METRO":
        from src.codon_opt.optimizers.metro_optimizer import MetropolisOptimizer

        MetropolisOptimizer(config)._optimize()
    elif config.args.codon_optimizer == "REMC":
        from src.codon_opt.optimizers.replica_exchange_mc import REMCOptimizer

        REMCOptimizer(config)._optimize()
    else:
        config.log.error(
            "Please select a valid optimizer. See python design.py -h for details."
        )
        raise NotImplementedError(
            "Please select a valid optimizer. See python design.py -h for details."
        )
