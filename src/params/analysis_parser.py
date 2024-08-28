import argparse
import os
import sys
import pickle
import datetime
from src.version.version import __version__
from src.logging.logging import Logging


class AnalysisParser(object):
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
        self._parse(args)
        log_obj = Logging()
        self.log = log_obj._create_log(
            self.__prog__, self.__version__, self.args.log_file_name
        )
        self._connect_to_db()
        self._validate()
        log_obj._log_args(self.log, arg_list=vars(self.args))
        self._query_details()

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
            help="Input (SQLite or Postgres) database from output of design.py",
        )
        self.parser.add_argument(
            "-at",
            "--analysis_type",
            default="fe_landscape",
            type=str,
            help="Specify type of analysis. Values: fe_landscape, fe_generation, fe_trajectory, codon_trajectory",
        )
        self.parser.add_argument(
            "-hv",
            "--hash_value",
            default=None,
            type=str,
            help="Hash value of an optimization. If none provided, the first optimization in the database will be used.",
        )
        self.parser.add_argument(
            "-l",
            "--log_file_name",
            default="quvax.log",
            type=str,
            help="Log file for recording certain output, warnings, and errors",
        )
        self.parser.add_argument(
            "-o",
            "--output",
            default="analysis_out.txt",
            type=str,
            help="Specify output file",
        )
        self.parser.add_argument(
            "-sd",
            "--random_seed",
            default=1,
            type=int,
            help="Random seed for sequence generation, optimization, and folding",
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

        if args is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(args)

    def _connect_to_db(self):
        if self.args.database_type == "sqlite":
            import sqlite3

            self.log.info("Connecting to database " + self.args.input)
            self.db = sqlite3.connect(self.args.input)
        elif self.args.database_type == "postgres":
            import psycopg2
            from configparser import ConfigParser

            # parse ini file
            parser = ConfigParser()
            parser.read(self.args.database_ini)
            ini_data = {_[0]: _[1] for _ in parser.items("postgresql")}
            self.db = psycopg2.connect(
                f"user={ini_data['user']} password={ini_data['password']} dbname={self.args.input}"
            )
        else:
            raise NotImplementedError(
                "Database type (-db) "
                + self.args.database
                + " not implemented. Options: sqlite, postgres."
            )
        self.db_cursor = self.db.cursor()

    def _validate(self):
        """
        Validate user input. TO DO: need to validate database structure?

        """

        return

    def _query_details(self):
        # query to get sim_detail columns info
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
        if self.args.hash_value:
            query = f"SELECT * FROM SIM_DETAILS WHERE hash_value = '{self.args.hash_value}';"
        else:
            query = f"SELECT * FROM SIM_DETAILS;"
        self.db_cursor.execute(query)
        sim_details = self.db_cursor.fetchall()
        if len(sim_details) == 0:
            self.log.error(
                "No data retrieved from database. Check your inputs or database structure."
            )
            raise ValueError(
                "No data retrieved from database. Check your inputs or database structure."
            )
        elif len(sim_details) > 1 and self.args.hash_value is None:
            self.log.info(
                "No hash value was specified and multiple optimizations are in the database. Using the first listed."
            )
        self.sim_details = {}  # dict to store details for later access
        self.log.info("Input Optimization Details:")
        for i in range(len(keys)):
            self.log.info(keys[i][0] + " : " + str(sim_details[0][i]))
            self.sim_details[keys[i][0]] = sim_details[0][i]
