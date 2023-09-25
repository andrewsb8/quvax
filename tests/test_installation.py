"""
Run from parent dir:
$ python -m tests.test_installation [-v INT]

"""

import unittest
import argparse

def _parse():
    '''
    Define command line arguments. Long options are used as variable names.
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbosity", type=int, default=2, help="Change amount of output in tests")

    return parser.parse_args()

if __name__ == '__main__':
    args = _parse()
    testsuite = unittest.TestLoader().discover('tests/test_cases/')
    unittest.TextTestRunner(verbosity=args.verbosity).run(testsuite)
