import argparse
from Bio.Seq import Seq
from Bio import SeqIO
import os
import sys
import logging
import datetime
import hashlib
from src.exceptions.exceptions import InvalidSequenceError


class DesignParser(object):
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
    n_trials : int
        Number of initial codon sequences to generate
    codon_optimizer : str
        Designation of outer loop optimizer
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
    checkpoint_interval : int
        Frequency to write checkpoint
    hash_value : str
        Hash used to identify optimizations within a database produced by design.py
    database_type : str
        String to choose database type to use for storing optimization data. Default: sqlite. Options: sqlite, postgres.
    database_ini : str
        Input file containing access information for postgres database.
    sequence_rejections : int,
        For use with MC optimizer only. Maximum number of rejections before a random sequence is proposed. Default: 3.
    num_sequence_changes : int,
        For use with MC optimizer only. Number of changes to propose for any given sequence. Default: 1.

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
        self.parser.add_argument(
            "-cc",
            "--convergence",
            default=0,
            type=int,
            help="Terminates optimization if new free energy minimum is not found within an integer number of generations.",
        )
        self.parser.add_argument(
            "-sr",
            "--sequence_rejections",
            default=3,
            type=int,
            help="For use with MC optimizer only. Maximum number of rejections before a random sequence is proposed. Default: 3.",
        )
        self.parser.add_argument(
            "-nc",
            "--num_sequence_changes",
            default=1,
            type=int,
            help="For use with MC optimizer only. Number of changes to propose for any given sequence. Default: 1.",
        )
        self.parser.add_argument(
            "-b",
            "--beta",
            default=1,
            type=float,
            help="For use with MC optimizer only. Value for 1/kT. Default: 1.",
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

        if self.args.codon_optimizer == "GA" and self.args.n_trials < 2:
            raise ValueError(
                "Population size (-n) for the genetic algorithm (-co GA) must be at least 2."
            )

        if self.args.codon_optimizer == "TFDE" and self.args.n_trials < 4:
            raise ValueError(
                "Population size (-n) for the TF differential evolutionary optimizer (-co TFDE) must be at least 4."
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

    def _connect_to_db(self, database):
        if self.args.database_type == "sqlite":
            import sqlite3

            db = sqlite3.connect(database)
        elif self.args.database_type == "postgres":
            import psycopg2
            from configparser import ConfigParser

            # parse ini file
            parser = ConfigParser()
            parser.read(self.args.database_ini)
            ini_data = {_[0]: _[1] for _ in parser.items("postgresql")}
            conn = psycopg2.connect(
                f"user={ini_data['user']} password={ini_data['password']} dbname=postgres"
            )
            cursor = conn.cursor()
            # try to create database, except will rollback
            try:
                conn.autocommit = True  # need to create database
                cursor.execute(f"CREATE DATABASE {database}")
                conn.commit()
            except:
                conn.rollback()
                cursor.close()
                self.log.info("Database exists in postgres client.")
            # connect to database
            db = psycopg2.connect(f"user={ini_data['user']} password={ini_data['password']} dbname={database}")
        else:
            raise NotImplementedError(
                "Database type (-db) "
                + self.args.database
                + " not implemented. Options: sqlite, postgres."
            )
        return db

    def _prepare_db(self):
        self.db = self._connect_to_db(self.args.output)
        self.db_cursor = self.db.cursor()
        if self.args.hash_value is None:
            self.args.hash_value = hashlib.shake_256(
                (str(datetime.datetime.now()) + self.seq).encode()
            ).hexdigest(5)
            self.log.info("Job Hash generated by hashlib: " + str(self.args.hash_value))
        try:
            if self.args.database_type == "sqlite":
                primary_key_type = "INTEGER"
            elif self.args.database_type == "postgres":
                primary_key_type = "SERIAL"
            self.db_cursor.execute(
                f"""CREATE TABLE SIM_DETAILS (sim_key {primary_key_type}
                PRIMARY KEY, protein_seq_file VARCHAR, protein_sequence VARCHAR,
                 target_sequence VARCHAR, generation_size INT,
                 codon_opt_iterations INT, optimizer VARCHAR(10),
                 random_seed INT, min_free_energy FLOAT,
                 target_min_free_energy FLOAT, rna_solver
                 VARCHAR(20), rna_folding_iterations INT, min_stem_len
                 INT, min_loop_len INT, species VARCHAR, coeff_max_bond INT,
                 coeff_stem_len INT, generations_sampled INT, state_file
                 VARCHAR, checkpoint_interval INT, convergence INT, hash_value VARCHAR,
                 sequence_rejections INT, num_sequence_changes INT, beta FLOAT);"""
            )
            self.db_cursor.execute(
                f"""CREATE TABLE OUTPUTS (index_key {primary_key_type}
                PRIMARY KEY, sim_key INT, population_key INT, generation INT,
                sequences VARCHAR, energies FLOAT, secondary_structure VARCHAR);"""
            )
            self.db_cursor.execute(
                f"""CREATE TABLE MFE_SEQUENCES (index_key {primary_key_type} PRIMARY KEY,
                sim_key INT, sequences VARCHAR, secondary_structure VARCHAR)"""
            )
            self.log.info("Created database tables in " + self.args.output + "\n\n")
        except:
            self.db.rollback()
            self.log.info("Adding data to existing tables within database.\n\n")
            self.db_cursor.execute(
                f"SELECT sim_key FROM SIM_DETAILS WHERE hash_value = '{self.args.hash_value}';"
            )
            if len(self.db_cursor.fetchall()) > 0:
                raise ValueError(
                    "Hash value already exists in database. Please specify another value."
                )
        self.db_cursor.execute(
            f"""INSERT INTO SIM_DETAILS (protein_seq_file, protein_sequence,
            target_sequence, generation_size, codon_opt_iterations, optimizer,
            random_seed, rna_solver, rna_folding_iterations, min_stem_len,
            min_loop_len, species, coeff_max_bond, coeff_stem_len, state_file,
            convergence, checkpoint_interval, hash_value, sequence_rejections,
            num_sequence_changes, beta) VALUES
            ('{self.args.input}', '{self.seq}', '{self.args.target}',
            '{self.args.n_trials}', '{self.args.codon_iterations}',
            '{self.args.codon_optimizer}', '{self.args.random_seed}',
            '{self.args.solver}', '{self.args.rna_iterations}',
            '{self.args.min_stem_len}', '{self.args.min_loop_len}',
            '{self.args.species}', '{self.args.coeff_max_bond}',
            '{self.args.coeff_stem_len}', '{self.args.state_file}',
            '{self.args.checkpoint_interval}', '{self.args.convergence}',
            '{self.args.hash_value}', '{self.args.sequence_rejections}',
            '{self.args.num_sequence_changes}', '{self.args.beta}');"""
        )
        self.db.commit()
        # retrieve the integer value of the key associated with the input protein sequence with associated hash value
        self.db_cursor.execute(
            f"SELECT sim_key FROM SIM_DETAILS WHERE hash_value = '{self.args.hash_value}';"
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

    def _load_db(self):
        """
        Function to load information from previous use of design.py when the
        --resume option is used. This function assumes information for only a
        single optimization is present in the input database

        """
        self.log.info("Loading info from database " + self.args.input)
        # need to pass self in below line because function is initiated from class method, no instance of class yet
        self.db = self._connect_to_db(self, self.args.input)
        self.db_cursor = self.db.cursor()

        if self.args.hash_value is not None:
            query = f"SELECT * FROM SIM_DETAILS WHERE hash_value = '{self.args.hash_value}';"
        else:
            query = f"SELECT * FROM SIM_DETAILS;"
        try:
            self.db_cursor.execute(query)
        except:
            self.log.error("There was an error retreiving data from the database.")
            raise ValueError("There was an error retreiving data from the database.")
        data = self.db_cursor.fetchall()

        if len(data) == 0:
            self.log.error(
                "No data retrieved from database. Check your inputs or database structure."
            )
            raise ValueError(
                "No data retrieved from database. Check your inputs or database structure."
            )
        elif len(data) > 1 and self.args.hash_value is None:
            self.log.info("No hash value was specified and multiple optimizations are in the database. Using the first listed.")

        print(vars(self))

        # manually assigning inputs from database
        self.sim_key = data[0][0]
        self.seq = data[0][2]
        if data[0][3] == "None":
            self.args.target = None
        else:
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
        self.args.convergence = data[0][19]
        self.args.checkpoint_interval = data[0][20]
        self.args.sequence_rejections = data[0][22]  # skip hash value
        self.args.num_sequence_changes = data[0][23]
        self.args.beta = data[0][24]

        # originally set the codon iterations to the original number set by user minus the number sampled in previous iterations
        self.args.codon_iterations = (
            self.args.codon_iterations - self.generations_sampled
        )
        # if original number of steps have been completed, and user extends the optimization
        if self.args.extend > 0:
            self.log.info(
                "Extending optimization by " + str(self.args.extend) + " steps"
            )
            self.args.codon_iterations += self.args.extend
            self.db_cursor.execute(
                f"UPDATE SIM_DETAILS SET codon_opt_iterations = '{self.args.codon_iterations + self.generations_sampled}' WHERE protein_sequence = '{self.seq}';"
            )
            self.db.commit()
        elif self.args.codon_iterations == 0 and self.args.extend == 0:
            raise ValueError(
                "Optimization complete. Use -e to extend the optimization if desired. See python design.py --resume -h for details."
            )
        elif self.args.extend < 0:
            raise ValueError("Value for -e cannot be less than zero.")

        # collect final generation of sequences from previous execution of design.py
        self.db_cursor.execute(
            f"SELECT sequences, energies from OUTPUTS WHERE sim_key = '{self.sim_key}' and generation = '{self.generations_sampled}';"
        )
        data = self.db_cursor.fetchall()
        self.initial_sequences = [data[i][0] for i in range(len(data))]
        self.energies = [data[j][1] for j in range(len(data))]
        self.log.info("Loaded info from database " + self.args.input)
