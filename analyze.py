from src.config.analysis_config import AnalysisConfig


if __name__ == "__main__":
    config = AnalysisConfig()

    if config.args.command == "fe_landscape":
        from src.analysis.analyses.fe_landscape import FreeEnergyLandscape

        FreeEnergyLandscape(config)
    elif config.args.command == "fe_generation":
        from src.analysis.analyses.fe_generation import FreeEnergyGeneration

        FreeEnergyGeneration(config)
    elif config.args.command == "fe_trajectory":
        from src.analysis.analyses.fe_trajectory import FreeEnergyTrajectory

        FreeEnergyTrajectory(config)
    elif config.args.command == "codon_trajectory":
        from src.analysis.analyses.codon_trajectory import CodonTrajectory

        CodonTrajectory(config)
    elif config.args.command == "compare_ct":
        from src.analysis.analyses.compare_ct import CompareCT

        CompareCT(config)
    elif config.args.command == "base_pair_types":
        from src.analysis.analyses.base_pair_types import BasePairTypes

        BasePairTypes(config)
    elif config.args.command == "compute_energy":
        from src.analysis.analyses.compute_energy import ComputeEnergy

        ComputeEnergy(config)
    elif config.args.command == "k_neighbor_energy":
        from src.analysis.analyses.k_neighbor_energy import KNeighborEnergySearch

        KNeighborEnergySearch(config)
    elif config.args.command == "base_pair_ranges":
        from src.analysis.analyses.base_pair_ranges import BasePairRanges

        BasePairRanges(config)
    elif config.args.command == "classify_stems":
        from src.analysis.analyses.classify_stems import ClassifyStems

        ClassifyStems(config)
    elif config.args.command == "contact_order":
        from src.analysis.analyses.contact_order import ContactOrder

        ContactOrder(config)
    elif config.args.command == "unfold":
        from src.analysis.analyses.unfold import Unfold

        Unfold(config)
    else:
        config.log.error(
            "Please select a valid analysis. See python analyze.py -h for details."
        )
        raise NotImplementedError(
            "Please select a valid analysis. See python analyze.py -h for details."
        )
