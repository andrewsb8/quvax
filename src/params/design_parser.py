import argparse
from Bio.Seq import Seq
from Bio import SeqIO
import os
import sys
import logging
import datetime
from src.exceptions.exceptions import InvalidSequenceError
import sqlite3


class DesignParser(object):
    """
    Parses command line inputs using argparse.

    Parser Options
    ----------
    input : str
        Input file name
    codon_iterations : int
        Iterations for codon optimizations (outer loop)
    rna_iterations : int
        Iterations for RNA folded energy calcuations (inner loop)
    n_trials : int
        Number of initial codon sequences to generate
    codon_optimizer : str
        Designation of outer loop optimizer
    min_stem_len : int
        Minimum number of stems required in RNA folding
    min_loop_len : int
        Minimum number of loops required in RNA folding
    solver : str
        Designation of solver for RNA folding
    coeff_max_bond : int
        Coefficient for maximizing the number of bonds in RNA folding
    coeff_stem_len : int
        Coefficient for energetically penalizing short stems in RNA folding
    log_file_name : str
        String for log file for writing program outputs, warnings, and errors
    species : str
        String to identify which species to generate codon frequencies
    output : str
        String to identify output sqlite database file
    random_seed : int
        Sets random seed for all optimizers and packages
    target : str
        Optional input to include target codon sequence

    """

    def __init__(self, args=None):
        self._parse(args)
        self._load_input()
        self._logging()
        self._validate()
        self._log_args()
        self._create_db()

    def _parse(self, args=None):
        """
        Define command line arguments. Long options are used as variable names.
        """
        self.__version__ = "QuVax v0.0.1"
        self.prog = "design.py"

        self.parser = argparse.ArgumentParser(
            prog=self.prog,
            description="QuVax: mRNA design guided by folding potential",
            epilog="Please report bugs to: https://github.com/andrewsb8/quvax/issues",
        )
        self.parser.add_argument(
            "--version", action="version", version=self.__version__
        )
        self.parser.add_argument(
            "-i", "--input", required=True, type=str, help="Input sequence"
        )
        self.parser.add_argument(
            "-c",
            "--codon_iterations",
            default=100,
            type=int,
            help="Number of codon optimization (outer loop) iterations",
        )
        self.parser.add_argument(
            "-r",
            "--rna_iterations",
            default=10000,
            type=int,
            help="Number of RNA folding (inner loop) iterations",
        )
        self.parser.add_argument(
            "-n",
            "--n_trials",
            default=10,
            type=int,
            help="Number of initial sequences generated",
        )
        self.parser.add_argument(
            "-co",
            "--codon_optimizer",
            default="TFDE",
            type=str,
            help="Options: Genetic Algorithm (GA), Tensorflow Differential Evolution (TFDE), Random Optimizer (RAND)",
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
            default="hybrid",
            type=str,
            help="Choice of solver for RNA folding. Options: hybrid",
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
            "-o",
            "--output",
            default="quvax.db",
            type=str,
            help="String to identify output sqlite database file. Default: quvax.db",
        )
        self.parser.add_argument(
            "-sd",
            "--random_seed",
            default=1,
            type=int,
            help="Random seed for sequence generation, optimization, and folding",
        )
        self.parser.add_argument(
            "-t",
            "--target",
            default=None,
            type=str,
            help="Optional input to include target codon sequence",
        )

        if args is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(args)

    def _load_input(self):
        self.seq = str(SeqIO.read(self.args.input, "fasta").seq)
        if self.args.target is not None:
            self.args.target = self.args.target.upper()

    def _logging(self):
        self.log = logging.getLogger(__name__)
        self.log.setLevel(logging.DEBUG)
        handler = logging.FileHandler(self.args.log_file_name, mode="w+")
        self.log.addHandler(handler)
        if os.path.isfile(self.args.log_file_name):
            logging.warning(
                "Log file "
                + self.args.log_file_name
                + " exists and will be overwritten."
            )
        self.log.info("Program Version : " + self.__version__)
        self.log.info("Execution Time : " + str(datetime.datetime.now()))
        self.log.info(
            "Command line: python " + self.prog + " " + " ".join(sys.argv[1:]) + "\n\n"
        )
        self.log.info("Warnings and Errors:\n")

    def _validate(self):
        """
        Validate user input.

        """

        self.seq = self.seq.upper()

        aas = "ACDEFGHIKLMNPQRSTVWY"

        if any(_ not in aas for _ in self.seq):
            raise InvalidSequenceError(
                "At least one character in the input sequence is invalid"
            )

        if set(self.seq).issubset(set("GCATU")):
            self.log.warning("Input protein sequence looks like an DNA sequence!")

        if self.args.target is not None:
            cs = "GCATU"

            if any(_ not in cs for _ in self.args.target):
                raise InvalidSequenceError("Target is not a codon sequence!")

            if len(self.args.target) != 3 * len(self.seq):
                raise ValueError(
                    "Target sequence is not the correct length to code for the input protein sequence!"
                )

            target_reshape = [
                [self.args.target[i : i + 3]]
                for i in range(0, len(self.args.target), 3)
            ]
            if target_reshape.count(["AUG"]) > self.seq.count("M"):
                raise ValueError(
                    "Your target sequence includes the start codon AUG but your input protein sequence does not contain amino acid M!"
                )

            if (
                any("UAA" in codon for codon in target_reshape)
                or any("UAG" in codon for codon in target_reshape)
                or any("UGA" in codon for codon in target_reshape)
            ):
                raise ValueError(
                    "Your target sequence includes stop codons UAG, UGA, or UAA!"
                )

        if self.args.codon_iterations < 1:
            raise ValueError(
                """
            --codon_iterations must be at least 1!

            """
            )

        if self.args.rna_iterations < 1:
            raise ValueError(
                """
            --rna_iterations must be at least 1!

            """
            )

        if self.args.n_trials < 1:
            raise ValueError(
                """
            --n_trials must be at least 1!

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

    def _log_args(self):
        self.log.info("\n\nList of Parameters:")
        self.log.info("Protein Sequence : " + self.seq)
        iterable_args = vars(self.args)
        for k in iterable_args:
            self.log.info(k + " : " + str(iterable_args[k]))
        self.log.info("\n\n")

    def _create_db(self):
        self.log.info("Creating database " + self.args.output)
        self.db = sqlite3.connect(self.args.output)
        self.db_cursor = self.db.cursor()
        # This will fail if a db already exists in this directory
        self.db_cursor.execute(
            f"CREATE TABLE SIM_DETAILS (sim_key INTEGER PRIMARY KEY, protein_sequence VARCHAR({len(self.seq)}), target_sequence VARCHAR({len(self.seq)*3}), generation_size INT UNSIGNED, number_generations INT UNSIGNED, optimizer VARCHAR(10), random_seed INT, min_free_energy FLOAT, target_min_free_energy FLOAT);"
        )
        # f strings do not work with INSERT statements
        self.db_cursor.execute(
            "INSERT INTO SIM_DETAILS (protein_sequence, target_sequence, generation_size, number_generations, optimizer, random_seed) VALUES (?, ?, ?, ?, ?, ?);",
            (
                self.seq,
                self.args.target,
                self.args.n_trials,
                self.args.codon_iterations,
                self.args.codon_optimizer,
                self.args.random_seed,
            ),
        )
        self.db_cursor.execute(
            f"CREATE TABLE OUTPUTS (index_key INTEGER PRIMARY KEY, sim_key INT UNSIGNED, population_key INT UNSIGNED, generation INT UNSIGNED, sequences VARCHAR({len(self.seq)*3}), energies FLOAT);"
        )
        self.db_cursor.execute(
            f"CREATE TABLE MFE_SEQUENCES (index_key INTEGER PRIMARY KEY, sequences VARCHAR({len(self.seq)*3}))"
        )
        self.db.commit()
        # retrieve the integer value of the key associated with the input protein sequence, there is no check for redundant sequences
        self.db_cursor.execute(
            f"SELECT sim_key FROM SIM_DETAILS WHERE protein_sequence = '{self.seq}';"
        )
        self.sim_key = self.db_cursor.fetchall()[0][0]
        self.log.info("Created database " + self.args.output + "\n\n")
