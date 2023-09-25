from src.qodon.initiate_sequences import GenerateInitialSequences
import argparse
from Bio.Seq import Seq
from Bio import SeqIO
import sys
import warnings

class Parser(object):
    '''
    Parses command line inputs using argparse.

    Parameters
    ----------
    input : string
        Input file name
    codon_iterations : int
        Iterations for codon optimizations (outer loop)
    rna_iterations : int
        Iterations for RNA folded energy calcuations (inner loop)
    n_trials : int
        Number of initial codon sequences to generate
    codon_optimizer : str
        Designation of outer loop optimizer
    min_stem_len : int
        Minimum number of stems required in RNA folding
    min_loop_len : int
        Minimum number of loops required in RNA folding
    solver : str
        Designation of solver for RNA folding
    coeff_max_bond : int
        Coefficient for maximizing the number of bonds in RNA folding
    coeff_stem_len : int
        Coefficient for energetically penalizing short stems in RNA folding

    '''
    def __init__(self):
        self._parse()
        self.seq = str(SeqIO.read(self.args.input,'fasta').seq)
        self._validate()

        codons = GenerateInitialSequences(self.seq, self.args.n_trials)
        self.code_map = codons.code_map
        self.initial_sequences = codons.initial_sequences
        del codons

    def _parse(self):
        '''
        Define command line arguments. Long options are used as variable names.
        '''

        self.parser = argparse.ArgumentParser(description='QuVax: mRNA design guided by folding potential',
                    epilog='Please report bugs to: https://github.com/andrewsb8/quvax/issues')
        self.parser.add_argument("-i", "--input", required=True, type=str, help="Input sequence")
        self.parser.add_argument("-c", "--codon_iterations", default=100, type=int, help="Number of codon optimization (outer loop) iterations")
        self.parser.add_argument("-r", "--rna_iterations", default=10000, type=int, help="Number of RNA folding (inner loop) iterations")
        self.parser.add_argument("-n", "--n_trials", default=10, type=int, help="Number of initial sequences generated")
        self.parser.add_argument("-co", "--codon_optimizer", default="TFDE", type=str, help="Options: Genetic Algorithm (GA), Tensorflow Differential Evolution (TFDE)")
        self.parser.add_argument("-ms", "--min_stem_len", default=3, type=int, help="Minimum length of a RNA stem")
        self.parser.add_argument("-ml", "--min_loop_len", default=3, type=int, help="Minimum length of a RNA loop")
        self.parser.add_argument("-s", "--solver", default='hybrid', type=str, help="Choice of solver for RNA folding. Options: hybrid")
        self.parser.add_argument("-cB", "--coeff_max_bond", default=1, type=int, help="Coefficient for term maximizing number of bonds")
        self.parser.add_argument("-cL", "--coeff_stem_len", default=10, type=int, help="Coefficient for term penalizing short stems")

        self.args = self.parser.parse_args()

    def _validate(self):
        '''
        Validate user input.

        '''

        #I can't think of a time where this will execute
        if not isinstance(self.seq,str):

            raise TypeError('''
            Input protein sequence must be a string! User provided
            input with type {}

            '''.format(type(self.seq)))

        self.seq = self.seq.upper()

        aas = "ACDEFGHIKLMNPQRSTVWY"

        if any(_ not in aas for _ in self.seq):
            print('Not a valid input sequence!')
            sys.exit(1)

        if set(self.seq).issubset(set('GCATU')):
            warnings.warn("Input protein sequence looks like an DNA sequence!")

        if self.args.codon_iterations < 1:
            raise ValueError('''
            --codon_iterations must be at least 1!

            ''')

        if self.args.rna_iterations < 1:
            raise ValueError('''
            --rna_iterations must be at least 1!

            ''')

        if self.args.n_trials < 1:
            raise ValueError('''
            --n_trials must be at least 1!

            ''')