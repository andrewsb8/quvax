import os
import sys
import logging
from src.params.analysis_parser import Parser


if __name__ == "__main__":
    config = Parser()

    if config.args == None:
        print("placeholder")
    else:
        config.log.error("Please select a valid optimizer. See python design.py -h for details.")
