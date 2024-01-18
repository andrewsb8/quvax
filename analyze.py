from src.params.analysis_parser import AnalysisParser


if __name__ == "__main__":
    config = AnalysisParser()

    if config.args.analysis_type == "fe-landscape":
        from src.analysis.analyses.fe_landscape import FreeEnergyLandscape
        FreeEnergyLandscape(config)
    else:
        config.log.error(
            "Please select a valid analysis. See python analyze.py -h for details."
        )
