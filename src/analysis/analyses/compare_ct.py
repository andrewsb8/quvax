from src.analysis.analysis import Analysis
import pandas as pd

class CompareCT(Analysis):
    """
    Compares input connectivity table with reference connectivity table. Outputs sensitivity, PPV,
    F1 score, specificity, number of bases correctly paired, and number of bases incorrectly paired.

    Parameters
    ----------
    config : AnalysisParser
        Object containing user inputs

    """

    def __init__(self, config):
        super().__init__(config)
        self._analyze()

    def _analyze(self):
        self.get_pairings(self.config.args.input, self.config.args.reference)
        truepos, trueneg, falsepos, falseneg = self.truthvalues(quvaxpairings, targetpairings)
        self.compare_tables(quvaxpairings, targetpairings)

    def ct_to_dataframe(self, ct_file):
        """
        Takes base connectivity data in the form of .ct files and converts it to a data frame.

        """
        df = pd.read_csv(ct_file, delim_whitespace=True, skiprows=1, header=None)
        df.columns = ["Index", "Nucleotide", "Previous", "Next", "Paired With", "Counter"]
        return df

    def get_pairings(self, quvax_ct_file, ref_ct_file):
        """
        Takes user-input .ct files and obtains base pairing data from it.

        """
        targetct = self.ct_to_dataframe(ref_ct_file)
        ref_pairings = targetct["Paired With"]
        quvaxct = self.ct_to_dataframe(quvax_ct_file)
        quvaxpairings = quvaxct["Paired With"]
        return

    def compare_tables(self, test_pairings, ref_pairings):
        """
        Inputs the data frame, outputs the sensitivity, PPV, F1 score,
        specificity, bases correctly paired, and bases with incorrect pairing.

        """
        print("Bases correctly paired: ", str(truepos))
        print("Bases with misidentified pairings (either missed or incorrect pairing): ", str(falsepos + falseneg))
        print("Sensitivity: ", str(self.sensitivity(truepos, falseneg)), " ; PPV: ", str(self.PPV(truepos, falsepos))," ; F1: ", str(self.F1(truepos, falsepos, falseneg)), " ; Specificity: ", str(self.specificity(trueneg, falsepos)))

    def sensitivity(self, truepos, falseneg):
        return truepos/(truepos+falseneg)

    def PPV(self, truepos, falsepos):
        return truepos/(truepos+falsepos)

    def F1(self, truepos, falsepos, falseneg):
        return 2*truepos/(2*truepos + falsepos + falseneg)

    def specificity(self, trueneg, falsepos):
        return trueneg/(trueneg+falsepos)

    def truthvalues(self, test_pairings, ref_pairings):
        """
        Gathers data from the tables to find how many true positives, true negatives,
        false positives, and false negatives there are given both data sets .

        """
        numtruepos = 0
        numtrueneg = 0
        numfalsepos = 0
        numfalseneg = 0
        for i in range(len(test_pairings)):
            if test_pairings[i] != 0 and test_pairings[i] == ref_pairings[i]:
                numtruepos += 1
            if test_pairings[i] == 0 and ref_pairings[i] == 0:
                numtrueneg += 1
            if test_pairings[i] != 0 and test_pairings[i] != ref_pairings[i]:
                numfalsepos += 1
            if test_pairings[i] == 0 and ref_pairings[i] != 0:
                numfalseneg += 1
        return numtruepos, numtrueneg, numfalsepos, numfalseneg
