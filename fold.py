from src.params.fold_parser import FoldParser


if __name__ == "__main__":
    config = FoldParser()

    if config.args.solver == "SA":
        from src.rna_folding.rna_folders.simulated_annealer import SimulatedAnnealer

        #after MC folder PR is merged, will not need to specify a sequence
        fold = SimulatedAnnealer(config.seq, config)
    else:
        config.log.error(
            "Please select a valid folder. See python fold.py -h for details."
        )

    #after MC folder PR is merged, can do the following
    #fold._fold(fold.config.seq)
    #print(fold.best_score)
    #eventually secondary structure
