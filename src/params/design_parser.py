import argparse
from Bio.Seq import Seq
from Bio import SeqIO
import os
import sys
import hashlib
import datetime
from src.exceptions.exceptions import InvalidSequenceError
from src.version.version import __version__
from src.logging.logging import Logging


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
        log_obj = Logging()
        self.log= log_obj._create_log(self.__prog__, self.__version__, self.args.log_file_name)
        self._load_input()
        self._validate()
        log_obj._log_args(self.log, arg_list=vars(self.args), protein_sequence=self.protein_sequence)
        self._prepare_db()

    @classmethod
    def _resume(cls, args=None):
        cls.__prog__ = "design.py"
        cls.__version__ = __version__
        cls._parse_resume(cls, args)
        log_obj = Logging()
        cls.log= log_obj._create_log(cls.__prog__, cls.__version__, cls.args.log_file_name)
        cls._load_db(cls)
        log_obj._log_args(cls.log, arg_list=vars(cls.args), protein_sequence=cls.protein_sequence)
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
            db = psycopg2.connect(
                f"user={ini_data['user']} password={ini_data['password']} dbname={database}"
            )
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
                (str(datetime.datetime.now()) + self.protein_sequence).encode()
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
                 target VARCHAR, population_size INT,
                 codon_iterations INT, codon_optimizer VARCHAR(10),
                 random_seed INT, min_free_energy FLOAT,
                 target_min_free_energy FLOAT, solver
                 VARCHAR(20), rna_iterations INT, min_stem_len
                 INT, min_loop_len INT, species VARCHAR, coeff_max_bond INT,
                 coeff_stem_len INT, generations_sampled INT, state_file
                 VARCHAR, checkpoint_interval INT, convergence INT, hash_value VARCHAR,
                 sequence_rejections INT, num_sequence_changes INT, beta FLOAT,
                 beta_max FLOAT, exchange_frequency INT, mutation_chance FLOAT,
                 crossover_probability FLOAT, convergence_count INT, span INT);"""
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
            target, population_size, codon_iterations, codon_optimizer,
            random_seed, solver, rna_iterations, min_stem_len,
            min_loop_len, species, coeff_max_bond, coeff_stem_len, state_file,
            convergence, checkpoint_interval, hash_value, sequence_rejections,
            num_sequence_changes, beta, beta_max, exchange_frequency,
            mutation_chance, crossover_probability, convergence_count, span) VALUES
            ('{self.args.input}', '{self.protein_sequence}', '{self.args.target}',
            '{self.args.population_size}', '{self.args.codon_iterations}',
            '{self.args.codon_optimizer}', '{self.args.random_seed}',
            '{self.args.solver}', '{self.args.rna_iterations}',
            '{self.args.min_stem_len}', '{self.args.min_loop_len}',
            '{self.args.species}', '{self.args.coeff_max_bond}',
            '{self.args.coeff_stem_len}', '{self.args.state_file}',
            '{self.args.convergence}', '{self.args.checkpoint_interval}',
            '{self.args.hash_value}', '{self.args.sequence_rejections}',
            '{self.args.num_sequence_changes}', '{self.args.beta}',
            '{self.args.beta_max}', '{self.args.exchange_frequency}',
            '{self.args.mutation_chance}', '{self.args.crossover_probability}', 0,
            '{self.args.span}');"""
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

        # get column/variable names
        if self.args.database_type == "sqlite":
            self.db_cursor.execute(
                f"SELECT name FROM pragma_table_info('SIM_DETAILS');"
            )
        elif self.args.database_type == "postgres":
            # table name must be lower case for postgres!
            self.db_cursor.execute(
                f"SELECT column_name FROM information_schema.columns WHERE table_name='sim_details' ORDER BY ordinal_position;"
            )
        keys = self.db_cursor.fetchall()
        keys = [_[0] for _ in keys]  # make list from list of tuples

        # get values of variables
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
            self.log.info(
                "No hash value was specified and multiple optimizations are in the database. Using the first listed."
            )

        data = list(data[0])  # make list from tuple

        # mapping values to member attributes
        # manually keeping track of values which are not considered cli arguments, so members of self not self.args
        not_args = [
            "sim_key",
            "protein_sequence",
            "min_free_energy",
            "generations_sampled",
            "target_min_free_energy",
        ]
        mapping = dict(zip(keys, data))
        for key, val in mapping.items():
            if val == "None":  # None values read as strings from db
                val = None
            if key in not_args:
                setattr(self, key, val)
            elif key not in not_args:
                setattr(self.args, key, val)
            else:
                raise ValueError(
                    f"({key}, {val}) undefined. Check your database structure. You could be using an older version of QuVax."
                )

        if self.args.convergence_count == self.args.convergence:
            self.log.info(
                "Optimization was converged in previous run. By resuming, the convergence counter will be reset to zero and the same convergence criteria will be used again."
            )
            self.args.convergence_count = 0

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
                f"UPDATE SIM_DETAILS SET codon_iterations = '{self.args.codon_iterations + self.generations_sampled}' WHERE sim_key = '{self.sim_key}';"
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
            f"SELECT sequences, energies, secondary_structure from OUTPUTS WHERE sim_key = '{self.sim_key}' and generation = '{self.generations_sampled}';"
        )
        data = self.db_cursor.fetchall()
        if len(data) != self.args.population_size:
            raise ValueError(
                "Data set retrieved from Outputs is larger than population size."
            )
        self.initial_sequences = [data[i][0] for i in range(len(data))]
        self.energies = [data[j][1] for j in range(len(data))]
        self.sec_structs = [data[k][2] for k in range(len(data))]
        self.log.info("Loaded info from database " + self.args.input)
