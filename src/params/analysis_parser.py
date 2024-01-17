import argparse
from Bio.Seq import Seq
from Bio import SeqIO
import os
import sys
import logging
import datetime


class AnalysisParser(object):
    '''
    Parses command line inputs using argparse.

    Parser Options
    ----------
    input : str
        Input file name
    analysis_type : str
        Specify which analysis to perform
    log_file_name : str
        String for log file for writing program outputs, warnings, and errors
    output : str
        String to identify output file containing information about optimization process
    random_seed : int
        Sets random seed for all optimizers and packages

    '''
    def __init__(self, args = None):
        self._parse(args)
        self._logging()
        self._validate()
        self._log_args()

    def _parse(self, args=None):
        '''
        Define command line arguments. Long options are used as variable names.
        '''
        self.__version__ = 'QuVax v0.0.1'
        self.prog = 'analysis.py'

        self.parser = argparse.ArgumentParser(prog=self.prog,
                                              description='QuVax: mRNA design guided by folding potential',
                                              epilog='Please report bugs to: https://github.com/andrewsb8/quvax/issues')
        self.parser.add_argument("--version", action='version', version=self.__version__)
        self.parser.add_argument("-i", "--input", required=True, type=str, help="Input file from output of design.py. (default .qu)")
        self.parser.add_argument("-at", "--analysis_type", default="fe-landscape", type=str, help="Specify type of analysis. Values: fe-landscape")
        self.parser.add_argument("-l", "--log_file_name", default="quvax.log", type=str, help="Log file for recording certain output, warnings, and errors")
        self.parser.add_argument("-o", "--output", default="analysis_out.txt", type=str, help="Specify output file.")
        self.parser.add_argument("-sd", "--random_seed", default=1, type=int, help="Random seed for sequence generation, optimization, and folding")

        if args == None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(args)

    def _logging(self):
        self.log = logging.getLogger(__name__)
        self.log.setLevel(logging.DEBUG)
        handler = logging.FileHandler(self.args.log_file_name, mode='w+')
        self.log.addHandler(handler)
        if os.path.isfile(self.args.log_file_name):
            logging.warning("Log file " + self.args.log_file_name + " exists and will be overwritten.")
        self.log.info("Program Version : " + self.__version__)
        self.log.info("Execution Time : " + str(datetime.datetime.now()))
        self.log.info("Command line: python " + self.prog + " " + ' '.join(sys.argv[1:]) + "\n\n")
        self.log.info("Warnings and Errors:\n")

    def _validate(self):
        '''
        Validate user input.

        '''

    def _log_args(self):
        self.log.info("\n\nList of Parameters:")
        self.log.info("Protein Sequence : " + self.seq)
        iterable_args = vars(self.args)
        for k in iterable_args:
            self.log.info(k + " : " + str(iterable_args[k]))
        self.log.info("\n\n")
