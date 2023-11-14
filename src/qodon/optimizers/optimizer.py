from abc import ABC, abstractmethod
from src.params.parser import Parser
from src.qodon.initiate_sequences import GenerateInitialSequences
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

        codons = GenerateInitialSequences(self.config.seq, self.config.args.n_trials)
        self.code_map = codons.code_map
        self.initial_sequences = codons.initial_sequences
        del codons

    @abstractmethod
    def _optimize(self):
        pass

    def _tf_fold(self, nseq):
        '''
        Compute Minimum Free Energy (MFE) of RNA fold.

        '''
        rna_ss = RNAFold(nseq, self.config)
        return rna_ss.best_score

    def _get_nc(self, res):
        '''
        Extract number of possible codons for each amino acid

        '''
        return len(self.code_map[res]['codons'])

    def _reverse_translate(self, members):
        '''
        Convert to nucleotide sequence from integer indices of code map

        '''

        get_seq = lambda se: ''.join([self.code_map[res]['codons'][se[i] % self._get_nc(res)] for i, res in enumerate(self.config.seq)])
        seqs = [get_seq(se) for se in members]
        return seqs

    def _verify_dna(self, sequence):
        '''
        Translate nucleotide sequence to make sure it matches input

        '''
        if self.config.seq != str(Seq(sequence).transcribe().translate()):
            raise ValueError(
                "Error: Codon sequence did not translate properly!")
