import os
import sys
import logging
from src.params.analysis_parser import AnalysisParser


if __name__ == "__main__":
    config = AnalysisParser()

    if config.args == None:
        print("placeholder")
    else:
        config.log.error("Please select a valid analysis. See python analyze.py -h for details.")
