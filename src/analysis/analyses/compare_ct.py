from src.analysis.analysis import Analysis
import pandas as pd

class CompareCT(Analysis):
    """
    Compares base pairs between test and reference connectivity table files. Outputs sensitivity, PPV,
    F1 score, specificity, bases correctly paired, and bases incorrectly paired.

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
        self._analyze()

    def _analyze(self):
        input_pairings, reference_pairings = self._get_pairings(self.config.args.input, self.config.args.reference)
        self._truth_values(input_pairings, reference_pairings)
        self._compare_cts(self.truepos, self.trueneg, self.falsepos, self.falseneg)

    def _ct_to_dataframe(self, ct_file):
        """
        Converts connectivity tables to dataframe

        """
        df = pd.read_csv(ct_file, delim_whitespace=True, skiprows=1, header=None)
        df.columns = ["Index", "Nucleotide", "Previous", "Next", "Paired With", "Counter"]
        return df

    def _get_pairings(self, input_ct_file, reference_ct_file):
        """
        Return base pair information from dataframe

        """
        reference_ct = self._ct_to_dataframe(reference_ct_file)
        ref_pairings = reference_ct["Paired With"]
        input_ct = self._ct_to_dataframe(input_ct_file)
        in_pairings = input_ct["Paired With"]
        return in_pairings, ref_pairings

    def _compare_cts(self, truepos, trueneg, falsepos, falseneg):
        """
        Inputs the data frame, outputs the sensitivity, PPV, F1 score,
        specificity, bases correctly paired, and bases with incorrect pairing.

        """
        print("Bases correctly paired: ", str(truepos))
        print("Bases with misidentified pairings (either missed or incorrect pairing): ", str(falsepos + falseneg))
        print("Sensitivity: ", self._sensitivity(truepos, falseneg), " ; PPV: ", self._pos_predict_val(truepos, falsepos)," ; F1: ", self._f1(truepos, falsepos, falseneg), " ; Specificity: ", self._specificity(trueneg, falsepos))
        if self._sensitivity(truepos, falseneg) == None:
            print("WARNING: There are no correct base pairing matches or missed pairs. Please check that your inputs are valid.")
        if self._pos_predict_val(truepos, falsepos) == None:
            print("WARNING: There are no base pairings in your input file.")
        if self._f1(truepos, falsepos, falseneg) == None:
            print("WARNING: Infeasible result, please be sure your inputs are valid files.")
        if self._specificity(trueneg, falsepos) == None:
            print("WARNING: Every base in the reference file is paired.")



    def _sensitivity(self, truepos, falseneg):
        try:
            return truepos/(truepos+falseneg)
        except ZeroDivisionError:
            return None

    def _pos_predict_val(self, truepos, falsepos):
        try:
            return truepos/(truepos+falsepos)
        except ZeroDivisionError:
            return None

    def _f1(self, truepos, falsepos, falseneg):
        try:
            return 2*truepos/(2*truepos + falsepos + falseneg)
        except ZeroDivisionError:
            return None

    def _specificity(self, trueneg, falsepos):
        try:
            return trueneg/(trueneg+falsepos)
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
            if test_pairings[i] == 0 and ref_pairings[i] == 0:
                self.trueneg += 1
            if test_pairings[i] != 0 and test_pairings[i] != ref_pairings[i]:
                self.falsepos += 1
            if test_pairings[i] == 0 and ref_pairings[i] != 0:
                self.falseneg += 1
        return
