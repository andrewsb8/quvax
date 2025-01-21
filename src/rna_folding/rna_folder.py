import itertools
from abc import ABC, abstractmethod
from src.config.config import Config
from src.rna_structure.structure import RNAStructure
from src.rna_structure.structure_io import StructureIO
from src.rna_structure.structure_convert import StructureConvert


class RNAFolder(ABC, RNAStructure, StructureIO, StructureConvert):
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
        self.interactions = self._get_all_interactions()
        self.twobody_penalty = 500000
        self.no_stem_penalty = 500000

    def _declare_stem_vars(self, sequence):
        self.nseq = sequence
        self.n = len(self.nseq)
        self.stems = []
        self.h = dict()
        self.J = dict()
        self._pairs = []

    def _fold_prep(self, sequence):
        self._declare_stem_vars(sequence)
        self.config.log.debug("Sequence Length: " + str(self.n))
        self._gen_stems()
        self.len_stem_list = len(self.stems)
        self.config.log.debug(
            "Finished generating stems. Number of possible stems: "
            + str(self.len_stem_list)
        )
        self._compute_h_and_J()
        self.config.log.debug("Finished generating Hamiltonian matrices.")

    def _gen_stems(self):
        """
        Generates a list of tuples with three entries:
        - first two elements are sequence indices of a base pair identifying the stem
        - third element is the length of the stem
        - Ex: self.stems = [(1, 13, 3), ()]
            self.stems[0] = (1, 13, 3) corresponds to stem with base pairs between
            nucleotides with index 1 and 13, 2 and 12, 3 and 11.

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

    def _compute_h_and_J(self):
        if self.len_stem_list == 0:
            # return "infinite" energy, no simulated annealing because no matrix
            # to build for the Hamiltonian
            self.best_score = self.no_stem_penalty
            return
        else:
            # Pull out stem lengths for simplicity
            stems = [_[2] for _ in self.stems]
            if self.config.args.target_stem_length != -1:
                mu = self.config.args.target_stem_length
            else:
                mu = max(stems)

        a = 10 # constant
        f = 0.1 # fraction of bits to select

        # Compute all local fields and couplings
        h = {
            ind: self.config.args.coeff_stem_len * (ki**2 - 2 * mu * ki + mu**2)
            - self.config.args.coeff_max_bond * ki**2 - 2*a*f*self.n
            for ind, ki in enumerate(stems)
        }
        J = {
            (ind1, ind2): -2 * self.config.args.coeff_max_bond * ki1 * ki2 + a
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
                    J[(i, j)] += self.config.args.pseudo_factor * abs(J[(i, j)])

        self.h = h
        self.J = J

    def _calc_score(self, idx):
        """
        Calculate the score for the current list of stems.

        TODO: This can be made cheaper with array broadcasting and smarter slicing

        Parameters
        ----------
        idx : list
            list of stem indices. returns 0 if idx is empty

        """

        idx.sort()
        score = sum([self.h[x] for x in idx])
        score = score + sum([self.J[x] for x in itertools.combinations(idx, 2)])

        return score

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
            self.connect_list = self._stems_to_connect_list(self.n, self.stems_used)
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
            self.connect_list = self._stems_to_connect_list(self.n, self.stems_used)
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
