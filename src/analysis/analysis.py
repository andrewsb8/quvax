import difflib
from abc import ABC, abstractmethod
from src.config.config import Config


class Analysis(ABC):
    """
    Parent class for all analyses. Contains methods used by multiple
    analyses.

    Parameters
    ----------
    config : AnalysisParser
        Object containing user inputs

    """

    def __init__(self, config: Config):
        self.config = config

    @abstractmethod
    def _analyze(self):
        pass

    def _print_output_2D(self, filename, values):
        """
        Prints user-specified output file with two columns which are not
        binned. What the two columns represent and what post-processing is
        necessary is then up to the user.

        Parameters
        ----------
        values = [ [column 1 numbers], [column 2 numbers] ]
        """
        out = open(filename, "w+")
        for i in range(len(values[0])):
            line = ""
            for j in range(len(values)):
                line += str(values[j][i]) + " "
            out.write(line + "\n")
        out.close()

    def _calc_codon_diff(self, ref_codon, codon):
        """
        Checks if codons in two sequences are different.

        """
        diff = difflib.context_diff(ref_codon, codon)
        for i, s in enumerate(diff):
            if s[0] == "+" or s[0] == "-":
                return 1
        return 0

    def _codon_diff_list(self, sequences, ref_sequence=None):
        if ref_sequence:
            return [
                sum(
                    [
                        self._calc_codon_diff(
                            ref_sequence[i * 3 : (i * 3) + 3],
                            sequences[j][i * 3 : (i * 3) + 3],
                        )
                        for i in range(int(len(sequences[0]) / 3))
                    ]
                )
                for j in range(len(sequences))
            ]
        else:
            return [
                sum(
                    [
                        self._calc_codon_diff(
                            sequences[j][i * 3 : (i * 3) + 3],
                            sequences[j + 1][i * 3 : (i * 3) + 3],
                        )
                        for i in range(int(len(sequences[0]) / 3))
                    ]
                )
                for j in range(len(sequences) - 1)
            ]

    def _calc_energy_diff(self, ref_energy, energy):
        return energy - ref_energy
