from src.params.design_parser import DesignParser


if __name__ == "__main__":
    config = DesignParser()

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
        config.log.error(
            "Please select a valid optimizer. See python design.py -h for details."
        )
