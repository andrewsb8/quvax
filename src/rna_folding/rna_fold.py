import numpy as np
import os, itertools
from src.params.parser import Parser
import warnings
warnings.filterwarnings("ignore")
# Visualization tool: http://rna.tbi.univie.ac.at/forna/


class RNAFold(object):
    '''
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
    failed : bool
        description
    best_combo : list
        description
    best_scroe : float
        description

    '''
    def __init__(self, nseq, config: Parser):
        self.config = config
        self.nseq = nseq #specify nseq here to avoid confusion with self.config.seq
        self.n = len(self.nseq)

        self.stems = []
        self.h = dict()
        self.J = dict()
        self._pairs = []
        self.interactions = [('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G'),
                             ('G', 'U'), ('U', 'G')]
        self.twobody_penalty = 500000
        self.pseudo_factor = 0.5
        self.no_stem_penalty = 500000
        self.execute()

        print(self.best_score)

    def execute(self):
        self._gen_stems()
        self._compute_h_and_J()

    def compute_dwave_sa(self):
        import neal
        sampler = neal.SimulatedAnnealingSampler()
        h2 = { (k,k) : v for k,v in self.h.items() }
        Q = self.J
        Q.update(h2)
        if len(self.stems) > 100:
            self.config.args.rna_iterations = self.config.args.rna_iterations*2
        sampleset = sampler.sample_qubo(Q, num_reads=10, num_sweeps=self.config.args.rna_iterations)
        self.stems_used = [_ for it,_ in enumerate(self.stems) if it in [k for k,v in sampleset.first.sample.items() if v==1]]
        self.best_score = sampleset.first.energy

    def compute_exact(self):
        self._find_best_combo()

    def _gen_stems(self):
        for i in range(self.n - 2 * self.config.args.min_stem_len - self.config.args.min_loop_len):
            for j in range(i + 2 * self.config.args.min_stem_len + self.config.args.min_loop_len - 1,
                           self.n):
                for k in range(self.n):
                    if i + k >= self.n: break
                    if (j-k) - (i+k) < self.config.args.min_loop_len: break
                    if (self.nseq[i + k], self.nseq[j - k]) in self.interactions:
                        if k >= self.config.args.min_stem_len - 1:  # len-1 because k starts from 0 (not 1)
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
            pair1 for pair1 in pairs1
            if not any(_ in pair1 for _ in np.array(pairs2).flatten())
        ]
        keep2 = [
            pair2 for pair2 in pairs2
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
        if stem_pair[first][0] <= stem_pair[second][0] <= stem_pair[first][
                1] and stem_pair[second][0] <= stem_pair[first][
                    1] <= stem_pair[second][1]:
            return True
        return False

    def _compute_h_and_J(self):
        #print(
        #    'Treatment of pseudoknots is questionable at best at the moment!')

        # Pull out stem lengths for simplicity
        stems = [_[2] for _ in self.stems]
        if len(stems) == 0:
            #return "infinite" energy, no simulated annealing because no matrix
            #to build for the Hamiltonian
            #print("HERE")
            self.best_score = self.no_stem_penalty
            return
        else:
            mu = max(stems)

        # Compute all local fields and couplings
        h = {
            ind: self.config.args.coeff_stem_len * (ki**2 - 2 * mu * ki + mu**2) - self.config.args.coeff_max_bond * ki**2
            for ind, ki in enumerate(stems)
        }
        J = {
            (ind1, ind2): -2 * self.config.args.coeff_max_bond * ki1 * ki2
             for ind1, ki1 in enumerate(stems)
             for ind2, ki2 in enumerate(stems) if ind2 > ind1
        }

        # Replace couplings with 'infinite' energies for clashes. Adjust couplings
        # in cases of pseudoknots.
        for i in range(len(self.stems)):
            for j in range(i + 1, len(self.stems)):

                # If there's overlap, add large penalty and continue
                overlap = self._detect_stem_overlap(self.stems[i],
                                                    self.stems[j])
                if overlap:
                    J[(i, j)] = self.twobody_penalty
                    continue

                # Check if pseudoknot
                is_pseudo = self._is_pseudo(self.stems[i], self.stems[j])

                if is_pseudo:
                    J[(i, j)] += self.pseudo_factor * abs(J[(i,j)])

        if len(self.stems) == 0:
            J = {(0, 1): 0}

        self.h = h
        self.J = J

        #self._compute_model()
        self.compute_dwave_sa()

    def _compute_model(self):
        arr = np.zeros((len(self.h),len(self.h)))
        for i in range(len(self.h)):
            for j in range(len(self.h)):
                if i == j:
                    arr[i][i] = self.h[i]
                elif i > j:
                    arr[j][i] = self.J[(j,i)]
                else:
                    arr[i][j] = self.J[(i,j)]
        self.model = arr

    def _find_best_combo(self):
        best_score = 1000000
        best_combo = []
        # Iterate through combinations of size sc (from 2 to N_stems)
        for sc in range(2, len(self.stems)):
            # Iterate through all combinations containing sc elements
            for combo in (itertools.combinations(list(range(len(self.stems))),
                                                 sc)):
                # Sum contributions from each individual stem
                onebody_cont = sum([self.h[_] for _ in combo])
                # Identify all possible 2-body interactions between these stems
                twobody_cont = sum(
                    [self.J[tb] for tb in itertools.combinations(combo, 2)])
                # Total score
                combo_score = onebody_cont + twobody_cont
                # Compare to 'best'
                if combo_score < best_score:
                    best_score = combo_score
                    best_combo = combo
        # If best_score > 0 then no non-overlapping combos were found
        if best_score > 0:
            # Find longest stem
            best_stem = max(self.stems, key=lambda x: x[2])
            # One-body term is stem length squared
            best_score = best_stem[2]**2
            # Combo list contains one stem only
            best_combo = [self.stems.index(best_stem)]
        self.best_score = best_score
        self.best_combo = best_combo
        self.stems_used = [self.stems[_stem] for _stem in self.best_combo]

    def _detailed_score(self):
        print('Stems used in SS:')
        self._stems_used = []
        for _stem in self.best_combo:
            print(_stem, self.stems[_stem])
            self._stems_used.append(self.stems[_stem])

        print('\nInteraction terms:')
        for ob in self.best_combo:
            print(ob, self.h[ob])
        for tb in itertools.combinations(self.best_combo, 2):
            print(tb, self.J[tb])

        print('\nTotal score:', self.best_score)

    def ss_output(self):
        if len(self.best_combo) == 0:
            print('Must compute exact solution before outputting SS')
            raise Exception

        stems_used = []
        for c in self.best_combo:
            stems_used.append(self.stems[c])

        lefts = []
        rights = []
        for i, j, k in stems_used:
            lefts.append([i + _ for _ in range(k)])
            rights.append([j - _ for _ in range(k)])
        lefts = np.array(lefts).flatten()
        rights = np.array(rights).flatten()

        # Output SS format
        ss_seq = []
        for pos in range(len(self.nseq)):
            if pos + 1 in lefts:
                ss_seq.append('(')
            elif pos + 1 in rights:
                ss_seq.append(')')
            else:
                ss_seq.append('.')

        self.ss_str = ''.join(ss_seq)
        print(self.nseq)
        print(self.ss_str)


if __name__ == "__main__":

    #seq = 'ACGCGGGUACUGCGAUAGUG'
    seq = 'ACGUGAAGGCUACGAUAGUGCCAG'
    ## BCRV1
    #seq = 'UAUAUACUAGGUUGGCAUUUUGAGCGCAUCUUACUCAAAUCCUAGUAUUUCCAUUAAUAUCUAAUGAUAUUAAUGAUGCCUCUUAAUAUAAGAGAUGC'
    rna_ss = RNAFold(seq)
    rna_ss.compute_exact()
    print('done')
    rna_ss._detailed_score()
    rna_ss.ss_output()
