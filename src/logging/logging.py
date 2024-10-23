import os
import sys
import logging
import datetime


class Log(object):
    def __init__(self):
        pass

    def _create_log(self, __prog__, __version__, log_file_name):
        log = logging.getLogger(__name__)
        log.setLevel(logging.DEBUG)
        if os.path.isfile(log_file_name):
            logging.warning(
                "Log file " + log_file_name + " exists and will be overwritten."
            )
        handler = logging.FileHandler(log_file_name, mode="w+")
        log.addHandler(handler)
        log.info("Program Version : " + __version__)
        log.info("Execution Time : " + str(datetime.datetime.now()))
        log.info(
            "Command line: python " + __prog__ + " " + " ".join(sys.argv[1:]) + "\n\n"
        )
        log.info("Warnings and Errors:\n")
        return log

    def _log_args(self, log, arg_list, protein_sequence=None):
        log.info("\n\nList of Parameters:")
        if protein_sequence is not None:
            log.info("Protein Sequence : " + protein_sequence)
        for k in arg_list:
            log.info(k + " : " + str(arg_list[k]))
        log.info("\n\n")
