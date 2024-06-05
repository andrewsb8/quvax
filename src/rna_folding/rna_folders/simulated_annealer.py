from src.rna_folding.rna_folder import RNAFolder


class SimulatedAnnealer(RNAFolder):
    """
    Find the mimimum energy folded structure of an RNA sequence using Simulated
    Annealing. The Hamiltonian and problem formulation allow for the folding
    energy determination to be run on quantum hardware. See Fox et al. PLoS,
    2022, https://doi.org/10.1371/journal.pcbi.1010032.


    Parameters
    ----------
    config : Parser
        Object containing user inputs
    stems_used : list
        List of stems in minimum energy folded structure
    best_score : float
        Lowest energy output from simulated annealer for RNA folding

    """

    def __init__(self, config):
        super().__init__(config)

    def _fold(self, sequence):
        self._fold_prep(sequence)
        if self.len_stem_list > 0:
            self._compute_dwave_sa()
            self._stems_to_dot_bracket(self.n, self.stems_used)

    def _compute_dwave_sa(self):
        import neal

        sampler = neal.SimulatedAnnealingSampler()
        h2 = {(k, k): v for k, v in self.h.items()}
        Q = self.J
        Q.update(h2)
        if self.len_stem_list > 100:
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
