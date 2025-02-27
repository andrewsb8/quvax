import argparse
import os
import sys
import pickle
import datetime
from src.config.config import Config
from src.version.version import __version__
from src.logging.logging import Log
from src.database.database import Database


class AnalysisConfig(Config):
    """
    Parses command line inputs using argparse.

    Parser Options
    ----------
    input : str
        Input file name
    analysis_type : str
        Specify which analysis to perform
    hash_value : int
        Hash used to identify optimizations within a database produced by design.py
    log_file_name : str
        String for log file for writing program outputs, warnings, and errors
    output : str
        String to identify output file containing information about optimization process
    random_seed : int
        Sets random seed for all optimizers and packages
    database_type : str
        String to choose database type to use for storing optimization data. Default: sqlite. Options: sqlite, postgres.
    database_ini : str
        Input file containing access information for postgres database.

    """

    def __init__(self, args=None):
        self.__prog__ = "analyze.py"
        self.__version__ = __version__
        self.args = self._parse(args)
        log_obj = Log()
        self.log = log_obj._create_log(
            self.__prog__, self.__version__, self.args.log_file_name
        )
        if hasattr(self.args, "database_type"):
            db_obj = Database()
            self.db, self.db_cursor = db_obj._connect_to_db(
                self.args.database_type,
                self.args.input,
                self.log,
                self.args.database_ini,
            )
            self.sim_details = db_obj._get_sim_details(self)
        self._validate()
        log_obj._log_args(self.log, arg_list=vars(self.args))

    def _parse(self, args=None):
        """
        Define command line arguments. Long options are used as variable names.
        """

        parser = argparse.ArgumentParser(
            prog=self.__prog__,
            description="Analysis tool with modules to analyze outputs of design.py and fold.py. There are tools to compare sequences, energies, and secondary structures. See descriptions of the available commands below.",
            epilog="Please report bugs to: https://github.com/andrewsb8/quvax/issues",
        )
        parser.add_argument("--version", action="version", version=self.__version__)

        # currently there are two types of analyses:
        # - analyses that take information from a database
        # - analyses that take information from a structure file or pair of structure files
        # going to create parsers with typical shared options for each group and then individual
        #   subparsers can have extra options which are unique to each analysis

        # database analysis shared options parser
        db_parser = argparse.ArgumentParser(add_help=False)
        db_parser.add_argument(
            "-i",
            "--input",
            required=True,
            type=str,
            help="Input (SQLite or Postgres) database from output of design.py",
        )
        db_parser.add_argument(
            "-hv",
            "--hash_value",
            default=None,
            type=str,
            help="Hash value of an optimization. If none provided, the first optimization in the database will be used.",
        )
        db_parser.add_argument(
            "-l",
            "--log_file_name",
            default="quvax.log",
            type=str,
            help="Log file for recording certain output, warnings, and errors",
        )
        db_parser.add_argument(
            "-o",
            "--output",
            default="analysis_out.txt",
            type=str,
            help="Specify output file",
        )
        db_parser.add_argument(
            "-db",
            "--database_type",
            default="sqlite",
            type=str,
            help="Option to choose database type. Default: sqlite. Options: sqlite, postgres.",
        )
        db_parser.add_argument(
            "-in",
            "--database_ini",
            default=None,
            type=str,
            help="database .ini file to connect to postgres database.",
        )

        # shared parser options for secondary structure analysis methods
        ss_parser = argparse.ArgumentParser(add_help=False)
        ss_parser.add_argument(
            "-i",
            "--input",
            required=True,
            type=str,
            help="Input secondary structure file (ex: connectivity table).",
        )
        ss_parser.add_argument(
            "-l",
            "--log_file_name",
            default="quvax.log",
            type=str,
            help="Log file for recording certain output, warnings, and errors",
        )

        # shared parser options for each analysis which involves computing energy
        energy_parser = argparse.ArgumentParser(add_help=False)
        energy_parser.add_argument(
            "-cB",
            "--coeff_max_bond",
            default=1,
            type=float,
            help="Coefficient for term maximizing number of bonds",
        )
        energy_parser.add_argument(
            "-cL",
            "--coeff_stem_len",
            default=10,
            type=float,
            help="Coefficient for term penalizing short stems",
        )
        energy_parser.add_argument(
            "-pf",
            "--pseudo_factor",
            default=0.5,
            type=float,
            help="Coefficient for term penalizing pseudknots. Default: 0.5",
        )
        energy_parser.add_argument(
            "-mu",
            "--target_stem_length",
            default=-1,
            type=float,
            help="Value defining target stem length for one-body penalty. Default: -1 (maximum possible stem)",
        )
        energy_parser.add_argument(
            "-A",
            "--equality_constraint_constant",
            default=0,
            type=float,
            help="Value defining constant for equality constraint. Default: 0",
        )
        energy_parser.add_argument(
            "-f",
            "--equality_constraint_fraction",
            default=0,
            type=float,
            help="Value defining fraction for equality constraint. Default: 0",
        )
        energy_parser.add_argument(
            "-ms",
            "--min_stem_len",
            type=int,
            default=3,
            help="Minimum stem length value",
        )
        energy_parser.add_argument(
            "-ml",
            "--min_loop_len",
            type=int,
            default=3,
            help="Minimum loop length value",
        )
        energy_parser.add_argument(
            "-sn",
            "--span",
            default=0,
            type=int,
            help="Option to specify maximum distance, in terms of relative sequence location, between base pairs that will be considered for stem formation. If < 1, no span will be used. Default: 0.",
        )

        # subparsers for each analysis which relies on a database
        subparsers = parser.add_subparsers(dest="command")
        parser_fe_landscape = subparsers.add_parser(
            "fe_landscape",
            parents=[db_parser],
            help="Calculate the relative 2D free energy landscape relative to the minimum energy sampled during design.",
        )
        parser_fe_trajectory = subparsers.add_parser(
            "fe_trajectory",
            parents=[db_parser],
            help="Print the evolution of the free energy over each generation for each population member.",
        )
        parser_fe_generation = subparsers.add_parser(
            "fe_generation",
            parents=[db_parser],
            help="Prints the relative free energy of a population member compared to the previous generation for each population member and for all generations.",
        )
        parser_codon_trajectory = subparsers.add_parser(
            "codon_trajectory",
            parents=[db_parser],
            help="Prints the number different codons in an RNA sequence relative to the previous generation for each population member and for all generations.",
        )
        parser_compare_ct = subparsers.add_parser(
            "compare_ct",
            parents=[ss_parser],
            help="Compares two connectivity tables by calculating metrics such as f1 score.",
        )
        parser_compare_ct.add_argument(
            "-r",
            "--reference",
            required=True,
            type=str,
            help="Reference secondary structure file (ex: connectivity table).",
        )
        parser_base_pair_ranges = subparsers.add_parser(
            "base_pair_ranges",
            parents=[ss_parser],
            help="Calculates the average, minimum, and maximum length, in number of bases, between base pairs for an input connectivity table.",
        )
        parser_classify_stems = subparsers.add_parser(
            "classify_stems",
            parents=[ss_parser],
            help="Calculates the average, minimum, and maximum stem length, sequence length, number of stems, number of stems in pseudoknots, and number of overlapping stems for an input connectivity table.",
        )
        parser_bpt = subparsers.add_parser(
            "base_pair_types",
            parents=[ss_parser],
            help="For given input structure file, calculates percent of bases in different types of secondary structures.",
        )
        parser_comp_energy = subparsers.add_parser(
            "compute_energy",
            parents=[ss_parser, energy_parser],
            help="For a given input sequence and structure (connectivity table or TODO dot-bracket), characterize energy landscape by changing k neighbors (stems).",
        )
        parser_kns = subparsers.add_parser(
            "k_neighbor_energy",
            parents=[ss_parser, energy_parser],
            help="For a given input sequence and structure (connectivity table or TODO dot-bracket), characterize energy landscape by changing k neighbors (stems).",
        )
        parser_kns.add_argument(
            "-k",
            "--neighbors",
            default=1,
            type=int,
            help="For a given input sequence and structure (connectivity table or TODO dot-bracket), calculate distribution of energies by changing k neighbors (stems)",
        )
        parser_contact_order = subparsers.add_parser(
            "contact_order",
            parents=[ss_parser],
            help="Calculate contact order of a given input structure",
        )
        parser_unfold = subparsers.add_parser(
            "unfold",
            parents=[ss_parser, energy_parser],
            help="Calculate incremental energies as stems are removed from structure",
        )
        parser_unfold.add_argument(
            "-sd",
            "--random_seed",
            default=1,
            type=int,
            help="Random seed which will change the order which stems are removed from the system.",
        )
        parser_stem_saturation = subparsers.add_parser(
            "stem_saturation",
            parents=[ss_parser, energy_parser],
            help="""For a given structure, determine how many stems can be added to structure. Note: not
            cB, cL, pf, and mu do not change the output of this tool.""",
        )

        if args is None:
            return parser.parse_args()
        else:
            return parser.parse_args(args)

    def _validate(self):
        """
        Validate user input. TO DO: need to validate database structure?

        """

        return
