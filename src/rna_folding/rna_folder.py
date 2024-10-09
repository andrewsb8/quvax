import numpy as np
from abc import ABC, abstractmethod
from src.config.config import Config
from src.rna_structure.structure_io import StructureIO
from src.rna_structure.structure_convert import StructureConvert


class RNAFolder(ABC, StructureIO, StructureConvert):
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
        Energetic penalty for overlapping stems
    no_stem_penalty : int
        Penalty for sequence with no possible stems for given min stem length/loop length
    pseudo_factor : float
        Penalty for pseudoknots

    """

    def __init__(self, config: Config):
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
                    if (j - k) - (i + k) < self.config.args.min_loop_len or (
                        self.config.args.span > 0
                        and (j - k) - (i + k) > self.config.args.span
                    ):
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

    def _post_process(self):
        """
        Method that will be called to log and write output files for
        the output of fold.py

        """
        self.config.log.info(
            "Folding energy of input codon sequence: " + str(self.best_score)
        )
        self.config.log.info("Folded secondary structure: " + str(self.dot_bracket))
        if self.config.args.output_type == "dot_bracket":
            self._write_dot_bracket(
                self.config.args.output, self.best_score, self.nseq, self.dot_bracket
            )
        elif self.config.args.output_type == "connect_table":
            self._write_connect_table(
                self.config.args.output,
                self.nseq,
                self.best_score,
                self.connect_list,
            )
        elif self.config.args.output_type == "all":
            self._write_dot_bracket(
                str(self.config.args.output + ".dot"),
                self.best_score,
                self.nseq,
                self.dot_bracket,
            )
            self._write_connect_table(
                str(self.config.args.output + ".ct"),
                self.nseq,
                self.best_score,
                self.connect_list,
            )
        else:
            raise ValueError(
                "Output type (-ot, --output_type) invalid. See python fold.py -h for details."
            )
