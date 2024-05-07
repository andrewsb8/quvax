from src.params.design_parser import DesignParser
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

    def __init__(self, config: DesignParser):
        super().__init__(config)

    def _fold(self, sequence):
        self._fold_prep(sequence)
        self._compute_dwave_sa()

    def _compute_dwave_sa(self):
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
