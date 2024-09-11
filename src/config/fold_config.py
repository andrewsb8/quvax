import argparse
from Bio.Seq import Seq
from Bio import SeqIO
import os
import sys
import datetime
from src.config.config import Config
from src.exceptions.exceptions import InvalidSequenceError
from src.logging.logging import Log
from src.version.version import __version__


class FoldConfig(Config):
    """
    Parses command line inputs using argparse.

    Parser Options
    ----------
    input : str
        Input codon sequence
    solver : str
        Designation of solver for RNA folding. Options: SA (Simulated Annealing), MC (Monte Carlo), ES (Exact Solver), CTSA (Co-transcriptional Simulated Annealing)
    min_stem_len : int
        Minimum number of stems required in RNA folding
    min_loop_len : int
        Minimum number of loops required in RNA folding
    coeff_max_bond : int
        Coefficient for maximizing the number of bonds in RNA folding
    coeff_stem_len : int
        Coefficient for energetically penalizing short stems in RNA folding
    span : int
        Option to specify maximum distance, in terms of relative sequence location, between base pairs that will be considered for stem formation. If < 1, no span will be used. Default: 0.
    log_file_name : str
        String for log file for writing program outputs, warnings, and errors
    species : str
        String to identify which species to generate codon frequencies
    random_seed : int
        Sets random seed for all optimizers and packages
    output : str
        File name to store folding energy, codon sequence, and secondary structure
    output_type : str
        Option to specify output type. Default: dot-bracket. Options: dot-bracket, connect_table, all.

    """

    def __init__(self, args=None):
        self.__prog__ = "fold.py"
        self.__version__ = __version__
        self.args = self._parse(args)
        log_obj = Log()
        self.log = log_obj._create_log(
            self.__prog__, self.__version__, self.args.log_file_name
        )
        self._load_input()
        self._validate()
        log_obj._log_args(self.log, arg_list=vars(self.args))

    def _parse(self, args=None):
        """
        Define command line arguments. Long options are used as variable names.
        """

        self.parser = argparse.ArgumentParser(
            prog=self.__prog__,
            description="Determine secondary structure of an RNA sequence",
            epilog="Please report bugs to: https://github.com/andrewsb8/quvax/issues",
        )
        self.parser.add_argument(
            "--version", action="version", version=self.__version__
        )
        self.parser.add_argument(
            "-i",
            "--input",
            required=True,
            type=str,
            help="Input codon sequence string",
        )
        self.parser.add_argument(
            "-r",
            "--rna_iterations",
            default=10000,
            type=int,
            help="Number of RNA folding (inner loop) iterations",
        )
        self.parser.add_argument(
            "-ms",
            "--min_stem_len",
            default=3,
            type=int,
            help="Minimum length of a RNA stem",
        )
        self.parser.add_argument(
            "-ml",
            "--min_loop_len",
            default=3,
            type=int,
            help="Minimum length of a RNA loop",
        )
        self.parser.add_argument(
            "-s",
            "--solver",
            default="SA",
            type=str,
            help="Choice of solver for RNA folding. Options: SA (Simulated Annealing), MC (Monte Carlo), ES (Exact Solver)",
        )
        self.parser.add_argument(
            "-cB",
            "--coeff_max_bond",
            default=1,
            type=int,
            help="Coefficient for term maximizing number of bonds",
        )
        self.parser.add_argument(
            "-cL",
            "--coeff_stem_len",
            default=10,
            type=int,
            help="Coefficient for term penalizing short stems",
        )
        self.parser.add_argument(
            "-sn",
            "--span",
            default=0,
            type=int,
            help="Option to specify maximum distance, in terms of relative sequence location, between base pairs that will be considered for stem formation. If < 1, no span will be used. Default: 0.",
        )
        self.parser.add_argument(
            "-l",
            "--log_file_name",
            default="quvax.log",
            type=str,
            help="Log file for recording certain output, warnings, and errors",
        )
        self.parser.add_argument(
            "-sp",
            "--species",
            default="h_sapiens_9606",
            type=str,
            help="Species type for generating codon tables and frequencies",
        )
        self.parser.add_argument(
            "-sd",
            "--random_seed",
            default=1,
            type=int,
            help="Random seed for sequence generation, optimization, and folding",
        )
        self.parser.add_argument(
            "-o",
            "--output",
            default="quvax.dot",
            type=str,
            help="String to name output dot-bracket file. Includes energy, sequence, and secondary structure. Default: quvax.dot (quvax.ct if --output_type is connect_table).",
        )
        self.parser.add_argument(
            "-ot",
            "--output_type",
            default="dot_bracket",
            type=str,
            help="Option to specify output type. Default: dot-bracket. Options: dot_bracket, connect_table, all.",
        )

        if args is None:
            return self.parser.parse_args()
        else:
            return self.parser.parse_args(args)

    def _load_input(self):
        # Folders take sequence string so no need to convert
        self.seq = self.args.input

    def _validate(self):
        """
        Validate user input.

        """

        self.seq = self.seq.upper()

        aas = "DEFHIKLMNPQRSTVWY"

        if any(_ in aas for _ in self.seq):
            raise InvalidSequenceError(
                "At least one character in the input sequence is invalid. Is this a DNA (T nucleotide) or protein sequence?"
            )

        # Could have a protein sequence of only glycine, cysteine, and alanine
        if set(self.seq).issubset(set("GCA")):
            self.log.warning("Input codon sequence looks like an protein sequence!")

        seq_reshape = [[self.seq[i : i + 3]] for i in range(0, len(self.seq), 3)]

        if self.args.span < 0:
            raise TypeError("Span (-sn) cannot be less than zero.")

        if self.args.span > len(self.seq):
            self.log.warning(
                "--span is longer than the codon sequence length. This is equivalent to span = 0. Check to make sure you used the correct value!"
            )

        if self.args.span != 0 and self.args.span < len(self.seq) * 0.3:
            self.log.warning(
                "--span is less than 30% of the sequence length. Low span value could prohibit secondary structure formation."
            )

        if self.args.rna_iterations < 1:
            raise ValueError(
                """
            --rna_iterations must be at least 1!

            """
            )

        if self.args.min_stem_len < 1:
            raise ValueError(
                """
            --min_stem_len must be at least 1!

            """
            )

        if self.args.min_loop_len < 1:
            raise ValueError(
                """
            --min_loop_len must be at least 1!

            """
            )
