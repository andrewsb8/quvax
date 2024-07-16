import random
from src.params.fold_parser import FoldParser


if __name__ == "__main__":
    config = FoldParser()
    random.seed(config.args.random_seed)

    if config.args.solver == "SA":
        from src.rna_folding.rna_folders.simulated_annealer import SimulatedAnnealer

        fold = SimulatedAnnealer(config)
    elif config.args.solver == "MC":
        from src.rna_folding.rna_folders.classical_mc import MC

        fold = MC(config)
    elif config.args.solver == "ES":
        from src.rna_folding.rna_folders.exact_solver import ExactSolver

        fold = ExactSolver(config)
    else:
        config.log.error(
            "Please select a valid folder. See python fold.py -h for details."
        )
        raise NotImplementedError(
            "Please select a valid folder. See python fold.py -h for details."
        )

    fold._fold(fold.config.seq)
    fold.config.log.info(
        "Folding energy of input codon sequence: " + str(fold.best_score)
    )
    fold.config.log.info("Folded secondary structure: " + str(fold.dot_bracket))

    output = open(fold.config.args.output, "w")
    output.write("> Folded energy: " + str(fold.best_score) + "\n")
    output.write(fold.config.seq + "\n")
    output.write(fold.dot_bracket + "\n")
    output.close()
