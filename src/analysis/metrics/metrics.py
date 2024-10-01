class Metrics(object):
    """
    Class for storing and calculating various metrics relevant for analyses

    """
    def __init__(self):
        self.truepos = 0
        self.trueneg = 0
        self.falsepos = 0
        self.falseneg = 0
        self.sensitivity = 0
        self.specificity = 0
        self.f1 = 0
        self.pos_predict_val = 0

    def _increment(self, num):
        num += 1

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
