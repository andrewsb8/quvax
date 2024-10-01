from src.analysis.analysis import Analysis
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
        self.truepos = 0
        self.trueneg = 0
        self.falsepos = 0
        self.falseneg = 0
        self.sensitivity = 0
        self.specificity = 0
        self.f1 = 0
        self.pos_predict_val = 0
        self._analyze()

    def _analyze(self):
        input_pairings = self._get_pairings(self.config.args.input)
        reference_pairings = self._get_pairings(self.config.args.reference)
        self._truth_values(input_pairings, reference_pairings)
        self._calculate_metrics(
            self.truepos, self.trueneg, self.falsepos, self.falseneg
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
        list = vars(self)
        for k in list:
            if list[k] == None:
                self.config.log.warning(
                    """Value for """ + k + """ is undefined. """
                    """Multiple truth values are zero which led to division by zero."""
                )
        self.config.log.info(", ".join([k for k in list if k != "config"]))
        vals = ", ".join([str(list[k]) for k in list if k != "config"])
        self.config.log.info(vals)
        # printing output values to stdout so it can be piped to file in bulk analysis
        print(vals)

    def _calculate_metrics(self, truepos, trueneg, falsepos, falseneg):
        """
        Calculates sensitivity, specificity, positive predictive value, F1 score

        """
        self.sensitivity = self._sensitivity(truepos, falseneg)
        self.specificity = self._specificity(trueneg, falsepos)
        self.f1 = self._f1(truepos, falsepos, falseneg)
        self.pos_predict_val = self._pos_predict_val(truepos, falsepos)

    def _sensitivity(self, truepos, falseneg):
        try:
            return truepos / (truepos + falseneg)
        except ZeroDivisionError:
            return None

    def _pos_predict_val(self, truepos, falsepos):
        try:
            return truepos / (truepos + falsepos)
        except ZeroDivisionError:
            return None

    def _f1(self, truepos, falsepos, falseneg):
        try:
            return 2 * truepos / (2 * truepos + falsepos + falseneg)
        except ZeroDivisionError:
            return None

    def _specificity(self, trueneg, falsepos):
        try:
            return trueneg / (trueneg + falsepos)
        except ZeroDivisionError:
            return None

    def _truth_values(self, test_pairings, ref_pairings):
        """
        Calculates true positives, true negatives,
        false positives, and false negatives from two sets of base pairs

        """
        for i in range(len(test_pairings)):
            if test_pairings[i] != 0 and test_pairings[i] == ref_pairings[i]:
                self.truepos += 1
            elif test_pairings[i] == 0 and ref_pairings[i] == 0:
                self.trueneg += 1
            elif test_pairings[i] != 0 and test_pairings[i] != ref_pairings[i]:
                self.falsepos += 1
            elif test_pairings[i] == 0 and ref_pairings[i] != 0:
                self.falseneg += 1
        return
