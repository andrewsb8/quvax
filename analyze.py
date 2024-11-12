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
    elif config.args.command == "sec_struct_types":
        from src.analysis.analyses.sec_struct_types import SecondaryStructureTypes

        SecondaryStructureTypes(config)
    elif config.args.command == "k_neighbor_energy":
        from src.analysis.analyses.k_neighbor_energy import kNeighborEnergySearch

        kNeighborEnergySearch(config)
    elif config.args.command == "pairing_ranges":
        from src.analysis.analyses.pairing_ranges import PairingRange

        PairingRange(config)
    else:
        config.log.error(
            "Please select a valid analysis. See python analyze.py -h for details."
        )
        raise NotImplementedError(
            "Please select a valid analysis. See python analyze.py -h for details."
        )
