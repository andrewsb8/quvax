import numpy as np
import os
import itertools
from src.params.design_parser import DesignParser
from src.rna_folding.rna_folder import RNAFolder


class QuantumSimAnnealer(RNAFolder):
    """
    Calculate the folded energy of a given codon sequence.

    Parameters
    ----------
    config : Parser
        Object containing user inputs
    nseq : list
        codon sequence for RNA folded energy calculation
    n : int
        length of codon sequence
    stems : list
        List of stems formed in folded RNA structure
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
    best_score : float
        Lowest energy output from simulated annealer for RNA folding

    """

    def __init__(self, nseq, config: DesignParser):
        self.config = config
        self.nseq = nseq  # specify nseq here to avoid confusion with self.config.seq
        self.n = len(self.nseq)

        self.stems = []
        self.h = dict()
        self.J = dict()
        self._pairs = []
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
        self.execute()

    def execute(self):
        self._gen_stems()
        self._compute_h_and_J()

    def compute_dwave_sa(self):
        import neal

        sampler = neal.SimulatedAnnealingSampler()
        h2 = {(k, k): v for k, v in self.h.items()}
        Q = self.J
        Q.update(h2)
        if len(self.stems) > 100:
            self.config.args.rna_iterations = self.config.args.rna_iterations * 2
        sampleset = sampler.sample_qubo(
            Q, num_reads=10, num_sweeps=self.config.args.rna_iterations
        )
        self.stems_used = [
            _
            for it, _ in enumerate(self.stems)
            if it in [k for k, v in sampleset.first.sample.items() if v == 1]
        ]
        self.best_score = sampleset.first.energy

    def _gen_stems(self):
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
        for i in range(len(self.stems)):
            for j in range(i + 1, len(self.stems)):
                # If there's overlap, add large penalty and continue
                overlap = self._detect_stem_overlap(self.stems[i], self.stems[j])
                if overlap:
                    J[(i, j)] = self.twobody_penalty
                    continue

                # Check if pseudoknot
                is_pseudo = self._is_pseudo(self.stems[i], self.stems[j])

                if is_pseudo:
                    J[(i, j)] += self.pseudo_factor * abs(J[(i, j)])

        if len(self.stems) == 0:
            J = {(0, 1): 0}

        self.h = h
        self.J = J

        self.compute_dwave_sa()
