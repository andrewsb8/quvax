import numpy as np
from random import uniform
import itertools, time, sys
from src.rna_folding.rna_folder import RNAFolder


class ExactSolver(RNAFolder):
    """
    Class to exactly minimize the rna folding hamiltonian for a given sequence

    """
    def __init__(self, config):
        super().__init__(config)

        ## MPI
        self.mpi_enabled = True
        self.comm = None
        self.size: int = None
        self.rank: int = None
        ## Check to see if we've got an MPI job
        self._init_mpi()

    def _fold(self, sequence):
        self._fold_prep(sequence)
        if 0 < self.len_stem_list < 30 or (self.len_stem_list > 30 and self.mpi_enabled):
            self._solve()
            self._stems_to_dot_bracket(self.n, self.stems_used)
        elif self.len_stem_list == 0:
            self._stems_to_dot_bracket(self.n, [])
        elif self.len_stem_list > 30 and not self.mpi_enabled:
            self.config.log.warning(f"{self.len_stem_list} stems detected. For systems with > 30 stems, it is highly recommended to use MPI. Otherwise, folding a single RNA sequence will take at least 1.5 hours and scale very poorly for longer sequences.")
        else:
            raise ValueError("Undefined behavior.")

    def _solve(self):
        self.best_score = 1e100
        self._dicts_to_np()
        ## Check if serial or parallel
        if self.size > 1:
            result = self._run_mpi()
        else:
            result = self._run_serial()

    def _init_mpi(self):
        """Checks to see if we're MPI enabled"""

        try:
            from mpi4py import MPI
        except Exception as e:
            self.comm = None
            self.rank = 0
            self.size = 1
            return
        else:
            self.comm = MPI.COMM_WORLD
            self.size = self.comm.Get_size()
            self.rank = self.comm.Get_rank()
            print(f"MPI Status: Rank: {self.rank} Size: {self.size}")
            sys.stdout.flush()

    def _run_serial(self):
        """Run solver serially"""

        ## Array of integers
        gen = itertools.count(self.rank, self.size)

        stop = 2**self.len_stem_list
        for num in gen:
            if num < stop:
                result, score = self._score(num)
                if score < self.best_score:
                    self.best_score = score
                    stem_indices = result.tolist()
            else:
                break

        self.stems_used = [
            self.stems[i] for i in range(len(stem_indices)) if stem_indices[i] == 1
        ]

    def _run_mpi(self):
        """Run solver with MPI"""

        ## Broadcast model/bits to other procs
        self.model = self.comm.bcast(self.model, root=0)

        ## Work on only the slice relevent to this rank
        gen = itertools.count(self.rank, self.size)
        stop = 2**self.len_stem_list
        for num in gen:
            if num < stop:
                result, score = self._score(num)
                if score < self.best_score:
                    best_score_ind = (score, result)
            else:
                break

        self.comm.barrier()
        best_score_all = self.comm.gather(best_score, root=0)
        if self.rank == 0:
            best_score_all = min(best_score_all, key=lambda x: x[1])
            self.best_score = best_score_all[0]
            stem_indices = best_score_all[1].tolist()
            self.stems_used = [
                self.stems[i] for i in range(len(stem_indices)) if stem_indices[i] == 1
            ]

    def _dicts_to_np(self):
        """Converts h/J dicts to 2D np array"""

        self.model = np.zeros((self.len_stem_list, self.len_stem_list), dtype="float32")
        for i in range(self.len_stem_list):
            for j in range(i, self.len_stem_list):
                if i == j:
                    self.model[i][i] = self.h[i]
                else:
                    self.model[i][j] = self.J[(i, j)]
        return self.model

    def _score(self, num):
        """Compute score for passed bit vector"""

        ## Convert integer to bit representation
        b = np.unpackbits(np.array([num], dtype=">i8").view(np.uint8))[-self.len_stem_list :]
        # b = np.array(list(np.binary_repr(num).zfill(self.len_stem_list))).astype(np.int8)
        b = b.reshape(-1, b.shape[0])
        result = np.einsum("ij,ik,jk->i", b, b, self.model)
        return b[result.argmin()], result.min()
