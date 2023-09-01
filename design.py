import sys
import argparse
from qodon.initiate_sequences import GenerateInitialSequences
from rna_folding.rna_fold import RNAFold
import tensorflow as tf
import numpy as np
import tensorflow_probability as tfp
from Bio.Seq import Seq
from Bio import SeqIO
import warnings


class QuDesign(object):

    def __init__(self):
        self._parse()
        self.seq = str(SeqIO.read(self.args.input,'fasta').seq)
        self._validate()

        codons = GenerateInitialSequences(self.seq, self.args.n_trials)
        self.code_map = codons.code_map
        self.initial_sequences = codons.initial_sequences
        del codons

        self._execute()

    def _parse(self):
        '''
        Define command line arguments. Long options are used as variable names.
        '''

        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("-i", "--input", required=True, type=str, help="Input sequence")
        self.parser.add_argument("-c", "--codon_iterations", default=100, type=int, help="Number of codon optimization (outer loop) iterations")
        self.parser.add_argument("-r", "--rna_iterations", default=10000, type=int, help="Number of RNA folding (inner loop) iterations")
        self.parser.add_argument("-n", "--n_trials", default=10, type=int, help="Number of initial sequences generated")

        self.args = self.parser.parse_args()

    def _validate(self):
        '''
        Validate user input.

        '''

        if not isinstance(self.seq,str):

            raise TypeError('''
            Input protein sequence must be a string! User provided
            input with type {}

            '''.format(type(self.seq)))

        self.seq = self.seq.upper()

        aas = 'ACDEFGHIKLMNPQRSTVWY'

        if any(_ not in aas for _ in self.seq):
            print('Not a valid input sequence!')
            exit()

        if set(self.seq).issubset(set('GCAU')):
            warnings.warn("Input protein sequence looks like an RNA sequence!")

        if set(self.seq).issubset(set('GCAT')):
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

    def _execute(self):
        '''
        Main execution. Run tensorflow optimizer for codon optimization. Objective
        function computes RNA structure with D-Wave's SA algorithm.

        '''

        self.initial_members = tf.convert_to_tensor(([_[1] for _ in self.initial_sequences]),np.float32)

        # Differential_weight: controls strength of mutations. We basically want to turn this off.
        # Crossover_prob: set this low. Need to think more about why this helps.
        optim_results = tfp.optimizer.differential_evolution_minimize(
            self._objective,
            initial_population=self.initial_members,
            max_iterations=self.args.codon_iterations,
            differential_weight=0.01,
            crossover_prob=0.1,
        )

        # Assign results as class attributes
        self.nseq = self._convert_to_nseqs(optim_results.final_population)[np.argmin(optim_results.final_objective_values)]
        self.mfe = np.min(optim_results.final_objective_values)

    def _objective(self,members):
        '''
        Objective function for TF to minimize

        NOTE: TF uses gradient descent to minimize continuous valued functions.
        The approach used here is not mathematically sound. It's a hack. But
        it gets the job done.

        '''

        # Map continuous valued tensor to RNA sequence
        n_seqs = self._convert_to_nseqs(members)

        # Use the imported scoring function to score all sequences.
        scores = [self._tf_fold(s) for s in n_seqs]

        # Return TF object
        return tf.cast(scores, np.float32)

    # Helper function to get number of possible codons for an amino acid
    def _get_nc(self,res):
        '''
        Extract number of possible codons for each amino acid

        '''
        return len(self.code_map[res]['codons'])

    def _tf_fold(self, nseq):
        '''
        Compute Minimum Free Energy (MFE) of RNA fold.

        '''
        rna_ss = RNAFold(nseq, min_stem_len=4, min_loop_len=4)
        results = rna_ss.compute_dwave_sa(sweeps=self.args.rna_iterations)
        return results.first.energy

    def _convert_to_nseqs(self, members):
        '''
        Continuous --> discrete transformation

        Doesn't make mathematical sense but it works.

        '''
        # This is a hack. TF deals with continuous valued functions. We need discrete and finite.
        # So let's cheat. Whatever values are assigned, make them ints and take the absolute value.
        members = np.absolute(np.array(members).astype(int))

        # Now we want to do something with the values. It's possible that some values exceed the
        # number of codons for the given position, so take the modulus. This is effectively a hashing
        # function. It's not mathematically rigorous, but it's good enough.
        # Finally, convert list of indices to the RNA sequence.
        get_seq = lambda se: ''.join([self.code_map[res]['codons'][se[i] % self._get_nc(res)] for i, res in enumerate(self.seq)])
        n_seqs = [get_seq(se) for se in members]
        return n_seqs


if __name__ == "__main__":
    exe = QuDesign()
    print(exe.mfe)
    print(exe.nseq)
