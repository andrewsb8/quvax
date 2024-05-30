from src.params.fold_parser import FoldParser


if __name__ == "__main__":
    config = FoldParser()

    if config.args.solver == "SA":
        from src.rna_folding.simulated_annealer import SimulatedAnnealer

        fold = SimulatedAnnealer(config)
    else:
        config.log.error(
            "Please select a valid folder. See python fold.py -h for details."
        )
