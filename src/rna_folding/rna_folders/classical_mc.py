from src.rna_folding.rna_folder import RNAFolder
import numpy as np
import random
import itertools
import copy
import math


class MC(RNAFolder):
    """
    Find the mimimum energy folded structure of an RNA sequence using Simulated
    Annealing. The Hamiltonian and problem formulation allow for the folding
    energy determination to be run on quantum hardware. See Fox et al. PLoS,
    2022, https://doi.org/10.1371/journal.pcbi.1010032.


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
            # self._log_mc_stats()

    def _add_pair(self):
        ## Grab a stem at random
        rand_idx = random.randint(0, self.len_stem_list - 1)

        ## Get the longest stem
        # rand_idx = self._get_largest_stem(rand_idx)

        if rand_idx in self.stem_idx:
            return

        stems = copy.copy(self.stem_idx)
        stems.append(rand_idx)

        ## Score the new set of interactions
        newscore = self._calc_score(stems)

        ## How'd we do?
        if newscore < self.score:
            self.stem_idx = stems
            self.score = newscore
            self.accept_add = self.accept_add + 1

        elif np.exp(-1 * (newscore - self.score) / self.T) > random.uniform(0.0, 1.0):
            self.stem_idx = stems
            self.score = newscore
            self.accept_add = self.accept_add + 1

        else:
            pass

    def _del_pair(self):
        # can't delete from set of zero
        if len(self.stem_idx) < 1:
            return

        ## delete a stem pair
        rand_idx = random.randint(0, len(self.stem_idx) - 1)
        stems = copy.copy(self.stem_idx)
        del stems[rand_idx]

        ## Score the new set of interactions
        newscore = self._calc_score(stems)

        ## How'd we do?
        if newscore < self.score:
            self.stem_idx = stems
            self.score = newscore
            self.accept_del = self.accept_del + 1

        elif np.exp(-1 * (newscore - self.score) / self.T) > random.uniform(0.0, 1.0):
            self.stem_idx = stems
            self.score = newscore
            self.accept_del = self.accept_del + 1

        else:
            pass

    def _swap_pair(self):
        # need at least two stems to swap
        if len(self.stem_idx) < 2:
            return

        rand_idx_del = random.randint(0, len(self.stem_idx) - 1)
        rand_idx_add = random.randint(0, self.len_stem_list - 1)

        stems = copy.copy(self.stem_idx)

        ## delete a stem pair
        del stems[rand_idx_del]

        ## Grab a stem at random
        # rand_idx_add = self._get_largest_stem(rand_idx_add)

        if rand_idx_add in self.stem_idx:
            return

        stems.append(rand_idx_add)

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

        else:
            pass  # why pass instead of return?

    # unused
    def _elongate_stem(self):
        """See if we can elongate a stem"""

        rand_idx = random.randint(0, len(self.stem_idx) - 1)
        stems = copy.copy(self.stem_idx)

        try:
            start, stop, length = self.stems[rand_idx]
            stems[rand_idx] = self.stems.index((start, stop, length + 1))
            newscore = self._calc_score(stems)
        except:
            ## Can't elongate stem based on pair list
            return

        ## How'd we do?
        if newscore < self.score:
            self.stem_idx = stems
            self.score = newscore
            # self.accept_swap = self.accept_swap + 1

        elif np.exp(-1 * (newscore - self.score) / self.T) > random.uniform(0.0, 1.0):
            self.stem_idx = stems
            self.score = newscore
            # self.accept_swap = self.accept_swap + 1

        else:
            pass

    # unused
    def _shorten_stem(self):
        """See if we can shorten a stem"""

        rand_idx = random.randint(1, len(self.stem_idx) - 1)
        stems = copy.copy(self.stem_idx)

        try:
            start, stop, length = self.stems[rand_idx]
            stems[rand_idx] = self.stems.index((start, stop, length - 1))
            newscore = self._calc_score(stems)
        except:
            ## Can't shorten stem based on pair list
            return

        ## How'd we do?
        if newscore < self.score:
            self.stem_idx = stems
            self.score = newscore
            # self.accept_swap = self.accept_swap + 1

        elif np.exp(-1 * (newscore - self.score) / self.T) > random.uniform(0.0, 1.0):
            self.stem_idx = stems
            self.score = newscore
            # self.accept_swap = self.accept_swap + 1

        else:
            pass

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
        selection of stem pair list. Sample size is set to math.floor(self.len_stem_list/5)
        to avoid small subset of stems for larger sequences. 5 was an arbitrary
        choice.

        TODO: Number of initial stems can probably be guessed by GC content.

        """

        self.stem_idx = random.sample(
            range(1, self.len_stem_list), math.floor(self.len_stem_list / 5)
        )
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

    def _do_mc(self, nsteps=100, T0=1.0):
        """Do some simple MC"""

        onethird: float = 1 / 3.0
        twothird: float = 2 / 3.0

        for inter in range(self.config.args.rna_iterations):
            ## Acceptance Ratio
            self.accept_add = 0
            self.accept_del = 0
            self.accept_swap = 0

            ## Initial Temperature
            self.T0 = T0

            ## Inital random guess for hairpins
            self._generate_init_ss_guess()

            for i in range(nsteps):
                ## Cool the system exponentially for now because it's easy
                self.T = self.T0 * np.exp(-i / nsteps)

                ## Choose a swap, insertion or deletion based on rando
                random_chance = random.uniform(0.0, 1.0)

                if random_chance <= onethird:
                    ## Attempt addition of stem pair
                    self._add_pair()
                elif onethird < random_chance <= twothird:
                    ## Attempt removal of stem pair
                    self._del_pair()
                else:
                    ## Attempt swap of stem pair from population
                    self._swap_pair()

            # self.best_score is returned to optimizer
            self.best_score = self.score

    # unused
    def _log_mc_stats(self):
        """
        May want to log these statistics but the output would be very long by
        default which could affect the usefulness of the log file with this folder.

        """

        print(f"Accept Ratio Add:  {self.accept_add  / float(nsteps)}")
        print(f"Accept Ratio Del:  {self.accept_del  / float(nsteps)}")
        print(f"Accept Ratio Swap: {self.accept_swap / float(nsteps)}")
        print(
            f"Accept Ratio: {(self.accept_add + self.accept_del + self.accept_swap) / float(nsteps)}"
        )
        print(f"Stems: {self.stem_idx}")
        print([self.stems[x] for x in self.stem_idx])

        for x in self.stem_idx:
            print(self.h[x])
        for x in itertools.combinations(self.stem_idx, 2):
            if x[0] < x[1]:
                print(self.J[x])
