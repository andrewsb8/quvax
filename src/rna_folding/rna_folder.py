import numpy as np
from abc import ABC, abstractmethod
from src.params.design_parser import DesignParser


class RNAFolder(ABC):
    """
    Parent class for all RNA folding classes.

    Parameters
    ----------
    nseq : list
        codon sequence for RNA folded energy calculation
    n : int
        length of codon sequence
    stems : list
        List of possible stems in folded RNA structure
    len_stem_list : int
        Length of list of stems
    h : dictionary
        Matrix for first term in folding Hamiltonian
    J : dictionary
        Matrix for second term in folding Hamiltonian
    interactions : list
        List of base pair interaction types
    pairs : list
        List of base pair interactions in the sequence
    twobody_penalty : int
        Energetic penalty
    pseudo_factor : float
        description

    """

    def __init__(self, config: DesignParser):
        self.config = config
        self.interactions = [
            ("A", "U"),
            ("U", "A"),
            ("G", "C"),
            ("C", "G"),
            ("G", "U"),
            ("U", "G"),
        ]
        self.twobody_penalty = 500000
        self.pseudo_factor = 0.5
        self.no_stem_penalty = 500000

    def _fold_prep(self, sequence):
        self.nseq = sequence
        self.n = len(self.nseq)
        self.stems = []
        self.h = dict()
        self.J = dict()
        self._pairs = []
        self._gen_stems()
        self.len_stem_list = len(self.stems)
        self._compute_h_and_J()

    def _gen_stems(self):
        """
        Generates a list of tuples with three entries:
        - first two elements are sequence indices of a base pair identifying the stem
        - third element is the length of the stem

        Can generate all base pair indices in a stem by calling _stem_to_pair_list
        Ex:
        - self.stems[0] -> (1, 13, 3)
        - _stem_to_pair_list(self.stems[0]) -> [(1, 13), (2, 12), (3, 11)]
        - Above shows the list of indices of three base pairs comprising the stem

        """

        for i in range(
            self.n - 2 * self.config.args.min_stem_len - self.config.args.min_loop_len
        ):
            for j in range(
                i
                + 2 * self.config.args.min_stem_len
                + self.config.args.min_loop_len
                - 1,
                self.n,
            ):
                for k in range(self.n):
                    if i + k >= self.n:
                        break
                    if (j - k) - (i + k) < self.config.args.min_loop_len:
                        break
                    if (self.nseq[i + k], self.nseq[j - k]) in self.interactions:
                        if (
                            k >= self.config.args.min_stem_len - 1
                        ):  # len-1 because k starts from 0 (not 1)
                            self._pairs.append((i + 1, j + 1, k + 1))
                    else:
                        break
        self.stems = self._pairs

    @staticmethod
    def _stem_to_pair_list(stem):
        pair_list = []
        for ci in range(stem[2]):
            pair_list.append((stem[0] + ci, stem[1] - ci))
        return pair_list

    def _detect_stem_overlap(self, stem1, stem2):
        pairs1 = self._stem_to_pair_list(stem1)
        pairs2 = self._stem_to_pair_list(stem2)

        keep1 = [
            pair1
            for pair1 in pairs1
            if not any(_ in pair1 for _ in np.array(pairs2).flatten())
        ]
        keep2 = [
            pair2
            for pair2 in pairs2
            if not any(_ in pair2 for _ in np.array(pairs1).flatten())
        ]

        if len(keep1) == len(pairs1) and len(keep2) == len(pairs2):
            # No overlap
            return False
        else:
            # Overlap
            return True

    def _is_pseudo(self, stem1, stem2):
        first = np.argmin([stem1[0], stem2[0]])
        second = np.argmax([stem1[0], stem2[0]])
        stem_pair = [stem1, stem2]
        if (
            stem_pair[first][0] <= stem_pair[second][0] <= stem_pair[first][1]
            and stem_pair[second][0] <= stem_pair[first][1] <= stem_pair[second][1]
        ):
            return True
        return False

    def _compute_h_and_J(self):
        # print(
        #    'Treatment of pseudoknots is questionable at best at the moment!')

        # Pull out stem lengths for simplicity
        stems = [_[2] for _ in self.stems]
        if len(stems) == 0:
            # return "infinite" energy, no simulated annealing because no matrix
            # to build for the Hamiltonian
            self.best_score = self.no_stem_penalty
            return
        else:
            mu = max(stems)

        # Compute all local fields and couplings
        h = {
            ind: self.config.args.coeff_stem_len * (ki**2 - 2 * mu * ki + mu**2)
            - self.config.args.coeff_max_bond * ki**2
            for ind, ki in enumerate(stems)
        }
        J = {
            (ind1, ind2): -2 * self.config.args.coeff_max_bond * ki1 * ki2
            for ind1, ki1 in enumerate(stems)
            for ind2, ki2 in enumerate(stems)
            if ind2 > ind1
        }

        # Replace couplings with 'infinite' energies for clashes. Adjust couplings
        # in cases of pseudoknots.
        for i in range(self.len_stem_list):
            for j in range(i + 1, self.len_stem_list):
                # If there's overlap, add large penalty and continue
                overlap = self._detect_stem_overlap(self.stems[i], self.stems[j])
                if overlap:
                    J[(i, j)] = self.twobody_penalty
                    continue

                # Check if pseudoknot
                is_pseudo = self._is_pseudo(self.stems[i], self.stems[j])

                if is_pseudo:
                    J[(i, j)] += self.pseudo_factor * abs(J[(i, j)])

        if self.len_stem_list == 0:
            J = {(0, 1): 0}

        self.h = h
        self.J = J

    def _stems_to_dot_bracket(self, sequence_len, stems):
        """
        Function to convert a list of stems in a sequence to a dot-bracket
        notation.

        """

        dot_bracket = ["." for i in range(sequence_len)]
        for stem in stems:
            stem_pair_list = self._stem_to_pair_list(stem)
            for i in range(len(stem_pair_list)):
                dot_bracket[stem_pair_list[i][0] - 1] = "("
                dot_bracket[stem_pair_list[i][1] - 1] = ")"

        # check for pseudoknots. pseudoknot cannot be detected if there is only one stem
        for i in range(len(stems)):
            for j in range(i + 1, len(stems)):
                if self._is_pseudo(stems[i], stems[j]):
                    stem_pair_list = self._stem_to_pair_list(stems[j])
                    for k in range(len(stem_pair_list)):
                        dot_bracket[stem_pair_list[k][0] - 1] = "["
                        dot_bracket[stem_pair_list[k][1] - 1] = "]"

        self.dot_bracket = "".join(dot_bracket)
