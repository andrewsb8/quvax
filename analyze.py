from src.params.analysis_parser import AnalysisParser


if __name__ == "__main__":
    config = AnalysisParser()

    if config.args.analysis_type == "fe_landscape":
        from src.analysis.analyses.fe_landscape import FreeEnergyLandscape

        FreeEnergyLandscape(config)
    elif config.args.analysis_type == "fe_generation":
        from src.analysis.analyses.fe_generation import FreeEnergyGeneration

        FreeEnergyGeneration(config)
    elif config.args.analysis_type == "trajectory":
        from src.analysis.analyses.trajectory import Trajectory

        Trajectory(config)
    elif config.args.analysis_type == "trajectory-pergeneration":
        from src.analysis.analyses.trajectory_per_generation import TrajectoryPerGeneration

        TrajectoryPerGeneration(config)
    else:
        config.log.error(
            "Please select a valid analysis. See python analyze.py -h for details."
        )
