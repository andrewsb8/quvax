from src.analysis.analysis import Analysis
from src.analysis.metrics.metrics import Metrics
import pandas as pd


class CompareCT(Analysis):
    """
    Compares base pairs between test and reference connectivity table files.
    Outputs truth values and sensitivity, specificity, PPV, F1 score.

    Parameters
    ----------
    config : AnalysisParser
        Object containing user inputs

    """

    def __init__(self, config):
        super().__init__(config)
        self.metrics = Metrics()
        self._analyze()

    def _analyze(self):
        input_pairings = self._get_pairings(self.config.args.input)
        reference_pairings = self._get_pairings(self.config.args.reference)
        self._truth_values(input_pairings, reference_pairings)
        self.metrics._calculate_metrics(
            self.metrics.truepos, self.metrics.trueneg, self.metrics.falsepos, self.metrics.falseneg
        )
        self._print_metrics()

    def _ct_to_dataframe(self, ct_file):
        """
        Converts connectivity table to dataframe

        """
        df = pd.read_csv(ct_file, delim_whitespace=True, skiprows=1, header=None)
        df.columns = [
            "Index",
            "Nucleotide",
            "Previous",
            "Next",
            "Paired With",
            "Counter",
        ]
        return df

    def _get_pairings(self, ct_file):
        """
        Return base pair information from dataframe

        """
        ct = self._ct_to_dataframe(ct_file)
        pairings = ct["Paired With"]
        return pairings

    def _print_metrics(self):
        list = vars(self.metrics)
        for k in list:
            if list[k] == None:
                self.config.log.warning(
                    """Value for """ + k + """ is undefined. """
                    """Multiple truth values are zero which led to division by zero."""
                )
        self.config.log.info(", ".join([k for k in list]))
        vals = ", ".join([str(list[k]) for k in list])
        self.config.log.info(vals)
        # printing output values to stdout so it can be piped to file in bulk analysis
        print(vals)

    def _truth_values(self, test_pairings, ref_pairings):
        """
        Calculates true positives, true negatives,
        false positives, and false negatives from two sets of base pairs

        """
        for i in range(len(test_pairings)):
            if test_pairings[i] != 0 and test_pairings[i] == ref_pairings[i]:
                self.metrics.truepos += 1
            elif test_pairings[i] == 0 and ref_pairings[i] == 0:
                self.metrics.trueneg += 1
            elif test_pairings[i] != 0 and test_pairings[i] != ref_pairings[i]:
                self.metrics.falsepos += 1
            elif test_pairings[i] == 0 and ref_pairings[i] != 0:
                self.metrics.falseneg += 1
        return
