from src.params.fold_parser import FoldParser


if __name__ == "__main__":
    config = FoldParser()

    if config.args.solver == "SA":
        from src.rna_folding.rna_folders.simulated_annealer import SimulatedAnnealer

        fold = SimulatedAnnealer(config)
    elif config.args.solver == "MC":
        from src.rna_folding.rna_folders.classical_mc import MC

        fold = MC(config)
    else:
        config.log.error(
            "Please select a valid folder. See python fold.py -h for details."
        )

    fold._fold(fold.config.seq)
    fold.config.log.info("Folding energy of input codon sequence: " + str(fold.best_score))

    output = open(fold.config.args.output, "w")
    output.write("> Folded energy: " + str(fold.best_score) + "\n")
    output.write(fold.config.seq)
    #output.write(secondary structure)
