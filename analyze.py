from src.params.analysis_parser import AnalysisParser


if __name__ == "__main__":
    config = AnalysisParser()

    if config.args.analysis_type == "fe_landscape":
        from src.analysis.analyses.fe_landscape import FreeEnergyLandscape

        FreeEnergyLandscape(config)
    elif config.args.analysis_type == "fe_generation":
        from src.analysis.analyses.fe_generation import FreeEnergyGeneration

        FreeEnergyGeneration(config)
    elif config.args.analysis_type == "fe_trajectory":
        from src.analysis.analyses.fe_trajectory import FETrajectory

        FETrajectory(config)
    elif config.args.analysis_type == "codon_trajectory":
        from src.analysis.analyses.codon_trajectory import CodonTrajectory

        CodonTrajectory(config)
    else:
        config.log.error(
            "Please select a valid analysis. See python analyze.py -h for details."
        )
