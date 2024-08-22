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

    fold._fold(fold.nseq)
    fold.config.log.info(
        "Folding energy of input codon sequence: " + str(fold.best_score)
    )
    fold.config.log.info("Folded secondary structure: " + str(fold.dot_bracket))
    if fold.config.args.output_type == "dot_bracket":
        fold._write_dot_bracket(
            fold.config.args.output, fold.best_score, fold.nseq, fold.dot_bracket
        )
    elif fold.config.args.output_type == "connect_table":
        fold._write_connect_table(
            fold.config.args.output,
            fold.nseq,
            fold.best_score,
            fold._stems_to_connect_list(fold.n, fold.stems_used),
        )
    elif fold.config.args.output_type == "all":
        fold._write_dot_bracket(
            str(fold.config.args.output + ".dot"),
            fold.best_score,
            fold.nseq,
            fold.dot_bracket,
        )
        fold._write_connect_table(
            str(fold.config.args.output + ".ct"),
            fold.nseq,
            fold.best_score,
            fold._stems_to_connect_list(fold.n, fold.stems_used),
        )
    else:
        raise ValueError(
            "Output type (-ot, --output_type) invalid. See python fold.py -h for details."
        )
