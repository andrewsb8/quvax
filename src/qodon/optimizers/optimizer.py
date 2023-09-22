from abc import ABC, abstractmethod
from src.params.parser import Parser
from src.rna_folding.rna_fold import RNAFold
from Bio.Seq import Seq

class Optimizer(ABC):
    """
    Parent class for all codon optimizer classes.

    Parameters
    ----------
    config : Parser
        Object containing user inputs

    """
    def __init__(self, config: Parser):
        self.config = config

    @abstractmethod
    def _optimize(self):
        pass

    def _tf_fold(self, nseq):
        '''
        Compute Minimum Free Energy (MFE) of RNA fold.

        '''
        rna_ss = RNAFold(nseq, self.config)
        results = rna_ss.compute_dwave_sa(sweeps=self.config.args.rna_iterations)
        return results.first.energy

    def _reverse_translate(self):
        '''
        Convert to nucleotide sequence

        '''
        self.final_codons = ''.join([
            self.config.code_map[res]['codons'][self.final_codons[i]]
            for i, res in enumerate(self.config.seq)
        ])

    def _verify_dna(self):
        '''
        Translate nucleotide sequence to make sure it matches input

        '''
        if self.config.seq != str(Seq(self.final_codons).transcribe().translate()):
            raise ValueError(
                "Error: Codon sequence did not translate properly!")
