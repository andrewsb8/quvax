import random
from src.config.fold_config import FoldConfig


if __name__ == "__main__":
    config = FoldConfig()
    random.seed(config.args.random_seed)

    if config.args.solver == "SA":
        from src.rna_folding.rna_folders.simulated_annealer import SimulatedAnnealer

        fold = SimulatedAnnealer(config)._fold(config.seq, post_process=True)
    elif config.args.solver == "CTSA":
        from src.rna_folding.rna_folders.cotranscript_sa import CoTranscriptSA

        fold = CoTranscriptSA(config)._fold(config.seq, post_process=True)
    elif config.args.solver == "MC":
        from src.rna_folding.rna_folders.classical_mc import MC

        fold = MC(config)._fold(config.seq, post_process=True)
    elif config.args.solver == "ES":
        from src.rna_folding.rna_folders.exact_solver import ExactSolver

        fold = ExactSolver(config)._fold(config.seq, post_process=True)
    else:
        config.log.error(
            "Please select a valid folder. See python fold.py -h for details."
        )
        raise NotImplementedError(
            "Please select a valid folder. See python fold.py -h for details."
        )
