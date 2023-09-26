from abc import ABC, abstractmethod
from src.params.parser import Parser
from src.rna_folding.rna_fold import RNAFold
from Bio.Seq import Seq
import numpy as np

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

    def _get_nc(self, res):
        '''
        Extract number of possible codons for each amino acid

        '''
        return len(self.config.code_map[res]['codons'])

    def _reverse_translate(self, members):
        '''
        Convert to nucleotide sequence from integer indices of code map

        '''

        get_seq = lambda se: ''.join([self.config.code_map[res]['codons'][se[i] % self._get_nc(res)] for i, res in enumerate(self.config.seq)])
        seqs = [get_seq(se) for se in members]
        return seqs

    def _verify_dna(self, sequence):
        '''
        Translate nucleotide sequence to make sure it matches input

        '''
        if self.config.seq != str(Seq(sequence).transcribe().translate()):
            raise ValueError(
                "Error: Codon sequence did not translate properly!")
