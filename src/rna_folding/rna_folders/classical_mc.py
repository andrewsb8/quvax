from src.rna_folding.rna_folder import RNAFolder
import numpy as np
import random
import itertools
import copy
import math


class MC(RNAFolder):
    """
    Find the mimimum energy folded structure of an RNA sequence using classical
    Monte Carlo. See Fox et al. PLoS, 2022, https://doi.org/10.1371/journal.pcbi.1010032.


    Parameters
    ----------
    config : Parser
        Object containing user inputs
    stems_idx : list
        List of stems in current folding state of RNA
    best_score : float
        Lowest energy output from simulated annealer for RNA folding

    """

    def __init__(self, config):
        super().__init__(config)

    def _fold(self, sequence):
        self._fold_prep(sequence)
        if self.len_stem_list > 0:
            self._do_mc()
            self._stems_to_dot_bracket(self.n, self.stems_used)
        else:
            self._stems_to_dot_bracket(self.n, [])

    def _add_pair(self):
        ## Grab a stem at random
        rand_idx = random.randint(0, self.len_stem_list - 1)

        if rand_idx in self.stem_idx:
            return self.stem_idx

        stems = copy.copy(self.stem_idx)
        stems.append(rand_idx)
        return stems

    def _del_pair(self):
        # can't delete from set of zero
        if len(self.stem_idx) < 1:
            return self.stem_idx

        ## delete a stem pair
        rand_idx = random.randint(0, len(self.stem_idx) - 1)
        stems = copy.copy(self.stem_idx)
        del stems[rand_idx]
        return stems

    def _swap_pair(self):
        # need at least two stems to swap
        if len(self.stem_idx) < 2:
            return self.stem_idx

        rand_idx_del = random.randint(0, len(self.stem_idx) - 1)
        rand_idx_add = random.randint(0, self.len_stem_list - 1)

        stems = copy.copy(self.stem_idx)

        ## delete a stem pair
        del stems[rand_idx_del]

        ## Grab a stem at random
        if rand_idx_add in self.stem_idx:
            return self.stem_idx

        stems.append(rand_idx_add)
        return stems

    def _get_largest_stem(self, stem_idx: int):
        """Given a stem index, get the largest stem from that start/stop group
        and return its corresponding index"""

        start, stop, length = self.stems[stem_idx]

        while (start, stop, length) in self.stems:
            length = length + 1

        return self.stems.index((start, stop, length - 1))

    def _generate_init_ss_guess(self):
        """
        Generate an initial plausable RNA starting point from random
        selection of stem pair list. Sample size is set to 3 so that


        TODO: Number of initial stems can probably be guessed by GC content.

        """

        self.stem_idx = random.sample(range(1, self.len_stem_list), 3)
        self.score = self._calc_score(self.stem_idx)

    def _calc_score(self, idx):
        """
        Calculate the score for the current list of stems

        TODO: This can be made cheaper with array broadcasting and smarter slicing

        """

        idx.sort()
        score = sum([self.h[x] for x in idx])
        score = score + sum([self.J[x] for x in itertools.combinations(idx, 2)])

        return score

    def _update_stems(self, stems):
        ## Score the new set of interactions
        newscore = self._calc_score(stems)

        ## How'd we do?
        if newscore < self.score:
            self.stem_idx = stems
            self.score = newscore
            self.accept_swap = self.accept_swap + 1

        elif np.exp(-1 * (newscore - self.score) / self.T) > random.uniform(0.0, 1.0):
            self.stem_idx = stems
            self.score = newscore
            self.accept_swap = self.accept_swap + 1

    def _do_mc(self, T0=1.0):
        """Do some simple MC"""

        onethird = 1 / 3.0
        twothird = 2 / 3.0

        ## Acceptance Ratio
        self.accept_add = 0
        self.accept_del = 0
        self.accept_swap = 0

        ## Initial Temperature
        self.T0 = T0

        ## Inital random guess for hairpins
        self._generate_init_ss_guess()

        for i in range(self.config.args.rna_iterations):
            ## Cool the system exponentially for now because it's easy
            self.T = self.T0 * np.exp(-i / self.config.args.rna_iterations)

            ## Choose a swap, insertion or deletion based on rando
            random_chance = random.uniform(0.0, 1.0)

            if random_chance <= onethird:
                ## Attempt addition of stem pair
                stems = self._add_pair()
            elif onethird < random_chance <= twothird:
                ## Attempt removal of stem pair
                stems = self._del_pair()
            else:
                ## Attempt swap of stem pair from population
                stems = self._swap_pair()

            if stems != self.stem_idx:
                self._update_stems(stems)

            # self.best_score is returned to optimizer
            if self.score < self.best_score:
                self.best_score = self.score
                self.stems_used = [self.stems[s] for s in self.stem_idx]
