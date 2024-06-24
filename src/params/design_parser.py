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
    state_file : str
        Output (or optional input file with --resume, see -h) to set random seed
        state
    hash_value : int
        Hash used to identify optimizations within a database produced by design.py

    """

    def __init__(self, args=None):
        self._parse(args)
        self._logging()
        self._load_input()
        self._validate()
        self._log_args()
        self._prepare_db()

    @classmethod
    def _resume(cls, args=None):
        cls._parse_resume(cls, args)
        cls._logging(cls)
        cls._load_db(cls)
        cls._log_args(cls)
        return cls

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
            default="SA",
            type=str,
            help="Choice of solver for RNA folding. Options: SA (Simulated Annealing), MC (Monte Carlo)",
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
            cs = "GCAU"

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

        if self.args.checkpoint_interval < 1:
            raise ValueError(
                """
            --checkpoint_interval must be at least 1!

            """
            )

        if self.args.checkpoint_interval > self.args.codon_iterations:
            self.log.warning(
                """
            Checkpoint interval is larger than the number of optimization steps!
            If you are running many steps, you might want to lower the
            checkpoint interval.

            """
            )

    def _log_args(self):
        self.log.info("\n\nList of Parameters:")
        self.log.info("Protein Sequence : " + self.seq)
        iterable_args = vars(self.args)
        for k in iterable_args:
            self.log.info(k + " : " + str(iterable_args[k]))
        self.log.info("\n\n")

    def _prepare_db(self):
        self.log.info("Creating database " + self.args.output)
        hash_value = hash(str(datetime.datetime.now()) + self.seq)
        self.log.info("Job Hash: " + str(hash_value))
        self.db = sqlite3.connect(self.args.output)
        self.db_cursor = self.db.cursor()
        try:
            self.db_cursor.execute(
                f"CREATE TABLE SIM_DETAILS (sim_key INTEGER PRIMARY KEY, protein_seq_file VARCHAR, protein_sequence VARCHAR, target_sequence VARCHAR, generation_size INT UNSIGNED, codon_opt_iterations INT UNSIGNED, optimizer VARCHAR(10), random_seed INT, min_free_energy FLOAT, target_min_free_energy FLOAT, rna_solver VARCHAR(20), rna_folding_iterations UNSIGNED INT, min_stem_len UNSIGNED INT, min_loop_len UNSIGNED INT, species VARCHAR, coeff_max_bond INT, coeff_stem_len INT, generations_sampled UNSIGNED INT, state_file VARCHAR, checkpoint_interval INT, hash_value INT);"
            )
            self.db_cursor.execute(
                f"CREATE TABLE OUTPUTS (index_key INTEGER PRIMARY KEY, sim_key INT UNSIGNED, population_key INT UNSIGNED, generation INT UNSIGNED, sequences VARCHAR, energies FLOAT, secondary_structure VARCHAR);"
            )
            self.db_cursor.execute(
                f"CREATE TABLE MFE_SEQUENCES (index_key INTEGER PRIMARY KEY, sim_key INT UNSIGNED, sequences VARCHAR({len(self.seq)*3}), secondary_structure VARCHAR)"
            )
            self.log.info("Created database " + self.args.output + "\n\n")
        except:
            self.log.info("Connected to existing database.\n\n")
        # f strings do not work with INSERT statements
        self.db_cursor.execute(
            "INSERT INTO SIM_DETAILS (protein_seq_file, protein_sequence, target_sequence, generation_size, codon_opt_iterations, optimizer, random_seed, rna_solver, rna_folding_iterations, min_stem_len, min_loop_len, species, coeff_max_bond, coeff_stem_len, state_file, checkpoint_interval, hash_value) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);",
            (
                self.args.input,
                self.seq,
                self.args.target,
                self.args.n_trials,
                self.args.codon_iterations,
                self.args.codon_optimizer,
                self.args.random_seed,
                self.args.solver,
                self.args.rna_iterations,
                self.args.min_stem_len,
                self.args.min_loop_len,
                self.args.species,
                self.args.coeff_max_bond,
                self.args.coeff_stem_len,
                self.args.state_file,
                self.args.checkpoint_interval,
                hash_value,
            ),
        )
        self.db.commit()
        # retrieve the integer value of the key associated with the input protein sequence, there is no check for redundant sequences
        self.db_cursor.execute(
            f"SELECT sim_key FROM SIM_DETAILS WHERE hash_value = '{hash_value}';"
        )
        self.sim_key = self.db_cursor.fetchall()[0][0]

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
            help="Input fasta-format protein sequence (or SQLite database with --resume)",
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
            type=int,
            help="Hash value of an optimization. If none provided, the first optimization in the database will be used.",
        )
        self.parser.add_argument(
            "--resume",
            action="store_true",
            help="Option to resume an optimization, -i needs to be a SQLite database file when using this flag and an input random state file is required for useful results",
        )

        if args is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(args)

    def _load_db(self):
        """
        Function to load information from previous use of design.py when the
        --resume option is used. This function assumes information for only a
        single optimization is present in the input database

        """
        self.log.info("Loading info from database " + self.args.input)
        self.db = sqlite3.connect(self.args.input)
        self.db_cursor = self.db.cursor()

        if self.args.hash_value is not None:
            query = f"SELECT * FROM SIM_DETAILS WHERE hash_value = '{self.args.hash_value}';"
        else:
            query = f"SELECT * FROM SIM_DETAILS;"
        try:
            self.db_cursor.execute(query)
        except:
            raise ValueError("Hash value not found in database.")
        data = self.db_cursor.fetchall()

        if len(data) == 0:
            raise ValueError("No data retrieved from database. Check your inputs.")

        # manually assigning inputs from database
        self.sim_key = data[0][0]
        self.seq = data[0][2]
        self.args.target = data[0][3]
        self.args.n_trials = data[0][4]
        self.args.codon_iterations = data[0][5]
        self.args.codon_optimizer = data[0][6]
        self.args.random_seed = data[0][7]
        self.mfe = data[0][8]
        self.target_folded_energy = data[0][9]
        self.args.solver = data[0][10]
        self.args.rna_iterations = data[0][11]
        self.args.min_stem_len = data[0][12]
        self.args.min_loop_len = data[0][13]
        self.args.species = data[0][14]
        self.args.coeff_max_bond = data[0][15]
        self.args.coeff_stem_len = data[0][16]
        self.generations_sampled = data[0][17]
        self.args.state_file = data[0][18]
        self.args.checkpoint_interval = data[0][19]
        if (
            self.args.hash_value is not None
        ):  # only overwrite it user did not provide a value
            self.args.hash_value = data[0][20]

        # originally set the codon iterations to the original number set by user minus the number sampled in previous iterations
        self.args.codon_iterations = (
            self.args.codon_iterations - self.generations_sampled
        )
        # if original number of steps have been completed, and user extends the optimization
        if self.args.codon_iterations == 0 and self.args.extend != 0:
            self.log.info(
                "Extending optimization by " + str(self.args.extend) + " steps"
            )
            self.args.codon_iterations += self.args.extend
            self.db_cursor.execute(
                "UPDATE SIM_DETAILS SET codon_opt_iterations = ? WHERE protein_sequence = ?;",
                (self.args.codon_iterations + self.generations_sampled, self.seq),
            )
            self.db.commit()
        elif self.args.codon_iterations == 0 and self.args.extend == 0:
            raise ValueError(
                "Optimization complete. Use -e to extend the optimization if desired. See python design.py --resume -h for details."
            )

        # collect final generation of sequences
        self.db_cursor.execute(
            f"SELECT sequences from OUTPUTS WHERE generation = {self.generations_sampled};"
        )
        sequences = self.db_cursor.fetchall()
        self.initial_sequences = [sequences[i][0] for i in range(len(sequences))]
        self.log.info("Loaded info from database " + self.args.input)
