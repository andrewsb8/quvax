import argparse
from Bio.Seq import Seq
from Bio import SeqIO
import os
import sys
from src.exceptions.exceptions import InvalidSequenceError
from src.version.version import __version__
from src.logging.logging import Logging
from src.database.database import Database


class DesignParser(Logging, Database):
    """
    Parses command line inputs using argparse.

    Parser Options
    ----------
    input : str
        Input file name
    codon_iterations : int
        Iterations for codon optimizations (outer loop)
    convergence : int
        Terminates optimization if new free energy minimum is not found within an integer number of generations
    rna_iterations : int
        Iterations for RNA folded energy calcuations (inner loop)
    population_size : int
        Number of codon sequences generated from input protein sequence to undergo optimization
    codon_optimizer : str
        Choice of outer loop optimizer. Ex: Metropolis Monte Carlo (METRO)
    min_stem_len : int
        Minimum number of stems required in RNA folding
    min_loop_len : int
        Minimum number of loops required in RNA folding
    solver : str
        Designation of solver for RNA folding. Options: SA (Simulated Annealing), MC (Monte Carlo), ES (Exact Solver)
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
    output : str
        String to identify output sqlite database file or postgres database
    random_seed : int
        Sets random seed for all optimizers and packages
    target : str
        Optional input to include a target codon sequence
    state_file : str
        Output (or optional input file with --resume, see -h) to set random seed
        state
    checkpoint_interval : int
        Interval of codon optimization steps to write checkpoint
    hash_value : str
        Hash used to identify optimizations within a database produced by design.py
    database_type : str
        String to choose database type to use for storing optimization data. Default: sqlite. Options: sqlite, postgres.
    database_ini : str
        Input file containing access information for postgres database.
    mutation_chance : float
        For use with GA optimizer only. Chance that a codon will randomly mutate in range [0, 1]. Default: 0.05.
    sequence_rejections : int,
        For use with METRO, REMC optimizers only. Maximum number of rejections before a random sequence is proposed. Default: 3.
    num_sequence_changes : int,
        For use with METRO, REMC optimizers only. Number of changes to propose for any given sequence. Default: 1.
    beta : float,
        For use with METRO, REMC optimizers only. Value for 1/kT to control "temperature" of optimization or acceptance probabilities. Lower beta (higher temperature) means changes are more likely to be accepted.
    beta_max : float,
        For use with REMC optimizer only. Max value of range for temperatures [beta, beta_max]. List of temperatures will be constructured with length equal to population size and interval (beta-max - beta)/population_size.
    exchange_frequency : int
        For use with REMC optimizer only. Frequency at which replica exchange attempts are made during codon optimization in generations.


    """

    def __init__(self, args=None):
        self.__prog__ = "design.py"
        self.__version__ = __version__
        self._parse(args)
        self.log = self._create_log(
            self.__prog__, self.__version__, self.args.log_file_name
        )
        self._load_input()
        self._validate()
        self._log_args(
            self.log, arg_list=vars(self.args), protein_sequence=self.protein_sequence
        )
        self.db, self.db_cursor = self._connect_to_db(self.args.database_type, self.args.output, self.log, self.args.database_ini, create=True)
        self._prepare_db(self)

    @classmethod
    def _resume(cls, args=None):
        cls.__prog__ = "design.py"
        cls.__version__ = __version__
        cls._parse_resume(cls, args)
        log_obj = Logging()
        cls.log = log_obj._create_log(
            cls.__prog__, cls.__version__, cls.args.log_file_name
        )
        db_obj = Database()
        cls.db, cls.db_cursor = db_obj._connect_to_db(cls.args.database_type, cls.args.input, cls.log, cls.args.database_ini)
        db_obj._load_db(cls)
        log_obj._log_args(
            cls.log, arg_list=vars(cls.args), protein_sequence=cls.protein_sequence
        )
        return cls

    def _parse(self, args=None):
        """
        Define command line arguments. Long options are used as variable names.
        """

        self.parser = argparse.ArgumentParser(
            prog=self.__prog__,
            description="QuVax: mRNA design guided by folding potential",
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
            help="Input fasta-format protein sequence (or SQLite database with --resume)",
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
            "-p",
            "--population_size",
            default=10,
            type=int,
            help="Number of codon sequences generated from input protein sequence to undergo optimization",
        )
        self.parser.add_argument(
            "-co",
            "--codon_optimizer",
            default="TFDE",
            type=str,
            help="Options: Genetic Algorithm (GA), Tensorflow Differential Evolution (TFDE), Random Optimizer (RAND), Metropolis Algorithm (METRO)",
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
            help="Output sqlite database file or postgres database. Default (sqlite): quvax.db",
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
        self.parser.add_argument(
            "-ci",
            "--checkpoint_interval",
            default=10,
            type=int,
            help="Frequency at which optimization details are updated: random state, min free energy, and generations sampled",
        )
        self.parser.add_argument(
            "--resume",
            action="store_true",
            help=argparse.SUPPRESS,
        )
        self.parser.add_argument(
            "-st",
            "--state_file",
            default="quvax.state",
            type=str,
            help="File to save (or load with --resume) the state of the pseudo random number generator",
        )
        self.parser.add_argument(
            "-cc",
            "--convergence",
            default=0,
            type=int,
            help="Terminates optimization if new free energy minimum is not found within an integer number of generations.",
        )
        self.parser.add_argument(
            "-sn",
            "--span",
            default=0,
            type=int,
            help="Option to specify maximum distance, in terms of relative sequence location, between base pairs that will be considered for stem formation. If < 1, no span will be used. Default: 0.",
        )
        self.parser.add_argument(
            "-cp",
            "--crossover_probability",
            default=0.10,
            type=float,
            help="For use with TFDE optimizer only. Probability of recombination for each codon in a sequence. Default: 0.10.",
        )
        self.parser.add_argument(
            "-mc",
            "--mutation_chance",
            default=0.05,
            type=float,
            help="For use with TFDE & GA optimizers only. Chance that a codon will randomly mutate in range [0, 1] ([0,2] for TFDE). Default: 0.05.",
        )
        self.parser.add_argument(
            "-sr",
            "--sequence_rejections",
            default=3,
            type=int,
            help="For use with METRO, REMC optimizers only. Maximum number of rejections before a random sequence is proposed. Default: 3.",
        )
        self.parser.add_argument(
            "-nc",
            "--num_sequence_changes",
            default=1,
            type=int,
            help="For use with METRO, REMC optimizers only. Number of changes to propose for any given sequence. Default: 1.",
        )
        self.parser.add_argument(
            "-b",
            "--beta",
            default=1,
            type=float,
            help="For use with METRO, REMC optimizers only. Value for 1/kT. Default: 1.",
        )
        self.parser.add_argument(
            "-bm",
            "--beta_max",
            default=10,
            type=float,
            help="For use with REMC optimizer only. Max value of range for temperatures [beta, beta_max]. List of temperatures will be constructured with length equal to population size and interval (beta-max - beta)/population_size. Default: 10.",
        )
        self.parser.add_argument(
            "-ef",
            "--exchange_frequency",
            default=10,
            type=int,
            help="For use with REMC optimizer only. Frequency at which replica exchange attempts are made during codon optimization in generations. Default: 10.",
        )
        self.parser.add_argument(
            "-db",
            "--database_type",
            default="sqlite",
            type=str,
            help="Option to choose database type. Default: sqlite. Options: sqlite, postgres.",
        )
        self.parser.add_argument(
            "-in",
            "--database_ini",
            default=None,
            type=str,
            help="database .ini file to connect to postgres database.",
        )
        self.parser.add_argument(
            "-hv",
            "--hash_value",
            default=None,
            type=str,
            help="Hash value to identify optimization in database. If none, one will be generated with hashlib. Default: None.",
        )

        if args is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(args)

    def _load_input(self):
        self.protein_sequence = str(SeqIO.read(self.args.input, "fasta").seq)
        if self.args.target is not None:
            self.args.target = self.args.target.upper()

    def _validate(self):
        """
        Validate user input.

        """

        self.protein_sequence = self.protein_sequence.upper()

        aas = "ACDEFGHIKLMNPQRSTVWY"

        if any(_ not in aas for _ in self.protein_sequence):
            raise InvalidSequenceError(
                "At least one character in the input sequence is invalid"
            )

        if set(self.protein_sequence).issubset(set("GCATU")):
            self.log.warning("Input protein sequence looks like an DNA sequence!")

        if self.args.span < 0:
            raise TypeError("Span (-sn) cannot be less than zero.")

        if self.args.span > len(self.protein_sequence) * 3:
            self.log.warning(
                "--span is longer than the codon sequence length. This is equivalent to span = 0. Check to make sure you used the correct value!"
            )

        if (
            self.args.span != 0
            and self.args.span < len(self.protein_sequence) * 3 * 0.3
        ):
            self.log.warning(
                "--span is less than 30% of the sequence length. Low span value could prohibit secondary structure formation."
            )

        if self.args.target is not None:
            cs = "GCAU"

            if any(_ not in cs for _ in self.args.target):
                raise InvalidSequenceError("Target is not a codon sequence!")

            if len(self.args.target) != 3 * len(self.protein_sequence):
                raise ValueError(
                    "Target sequence is not the correct length to code for the input protein sequence!"
                )

            target_reshape = [
                [self.args.target[i : i + 3]]
                for i in range(0, len(self.args.target), 3)
            ]
            if target_reshape.count(["AUG"]) != self.protein_sequence.count("M"):
                raise ValueError(
                    "Your target RNA sequence and protein sequence contain an unequal number of AUG codons and M amino acid residues! Did you leave the start codon in your target sequence?"
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

        if self.args.population_size < 1:
            raise ValueError(
                """
            --population_size must be at least 1!

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

        if self.args.checkpoint_interval < 1:
            raise ValueError(
                """
            --checkpoint_interval must be at least 1!

            """
            )

        if self.args.codon_optimizer == "GA" and self.args.population_size < 2:
            raise ValueError(
                "Population size (-p) for the genetic algorithm (-co GA) must be at least 2."
            )

        if self.args.codon_optimizer == "TFDE" and self.args.population_size < 4:
            raise ValueError(
                "Population size (-p) for the TF differential evolutionary optimizer (-co TFDE) must be at least 4."
            )

        if (
            self.args.codon_optimizer == "METRO" or self.args.codon_optimizer == "REMC"
        ) and (
            self.args.num_sequence_changes < 1
            or self.args.num_sequence_changes > len(self.protein_sequence)
        ):
            raise ValueError(
                "The number of changes proposed to any codon sequence (-nc) must be > 0 and <= the length of the input protein sequence."
            )

        if self.args.checkpoint_interval > self.args.codon_iterations:
            self.log.warning(
                """
            Checkpoint interval is larger than the number of optimization steps!
            If you are running many steps, you might want to lower the
            checkpoint interval.

            """
            )

        if not 0 <= self.args.mutation_chance <= 1:
            raise ValueError("mutation_chance (-mc) must be in the range [0, 1].")

        if self.args.beta <= 0:
            raise ValueError("beta (-b) cannot be negative.")

        if self.args.beta > self.args.beta_max:
            raise ValueError("beta_max (-bm) must be larger than beta (-b).")

    def _parse_resume(self, args=None):
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
            "-i",
            "--input",
            required=True,
            type=str,
            help="Database with information to resume optimization",
        )
        self.parser.add_argument(
            "-e",
            "--extend",
            default=0,
            type=int,
            help="Option to extend optimization by integer number of steps",
        )
        self.parser.add_argument(
            "-l",
            "--log_file_name",
            default="quvax.log",
            type=str,
            help="Log file for recording certain output, warnings, and errors",
        )
        self.parser.add_argument(
            "-st",
            "--state_file",
            default="quvax.state",
            type=str,
            help="File to save (or load with --resume) the state of the pseudo random number generator",
        )
        self.parser.add_argument(
            "-hv",
            "--hash_value",
            default=None,
            type=str,
            help="Hash value of an optimization. If none provided, the first optimization in the database will be used.",
        )
        self.parser.add_argument(
            "--resume",
            action="store_true",
            help="Option to resume an optimization, -i needs to be a SQLite database file when using this flag and an input random state file is required for useful results",
        )
        self.parser.add_argument(
            "-db",
            "--database_type",
            default="sqlite",
            type=str,
            help="Option to choose database type to retrieve optimization from. Default: sqlite. Options: sqlite, postgres.",
        )
        self.parser.add_argument(
            "-in",
            "--database_ini",
            default=None,
            type=str,
            help="database .ini file to connect to postgres database.",
        )

        if args is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(args)
