from qodon.optimizers.optimizer import Optimizer
import numpy as np
import tensorflow as tf
import tensorflow_probability as tfp
from rna_folding.rna_fold import RNAFold

class TfDiffEv(Optimizer):
    '''
    Tensorflow Differential Evolution optimizer for codon optimization.

    '''
    def __init__(self, config):
        super().__init__(config)
        self._optimize()

    def _optimize(self):
        '''
        Main execution. Run tensorflow optimizer. Objective
        function computes RNA structure with D-Wave's SA algorithm.

        '''

        self.initial_members = tf.convert_to_tensor(([_[1] for _ in self.config.initial_sequences]),np.float32)

        # Differential_weight: controls strength of mutations. We basically want to turn this off.
        # Crossover_prob: set this low. Need to think more about why this helps.
        optim_results = tfp.optimizer.differential_evolution_minimize(
            self._objective,
            initial_population=self.initial_members,
            max_iterations=self.config.args.codon_iterations,
            differential_weight=0.01,
            crossover_prob=0.1,
        )

        # Assign results as class attributes
        self.nseq = self._convert_to_nseqs(optim_results.final_population)[np.argmin(optim_results.final_objective_values)]
        self.mfe = np.min(optim_results.final_objective_values)

        print(self.mfe)
        print(self.nseq)

    def _objective(self, members): #without access to "parent" class variable values, need to copy QuDesign variables
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

    def _get_nc(self, res):
        '''
        Extract number of possible codons for each amino acid

        '''
        return len(self.config.code_map[res]['codons']) #HERE

    def _tf_fold(self, nseq):
        '''
        Compute Minimum Free Energy (MFE) of RNA fold.

        '''
        rna_ss = RNAFold(nseq, self.config)
        results = rna_ss.compute_dwave_sa(sweeps=self.config.args.rna_iterations) #HERE
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
        get_seq = lambda se: ''.join([self.config.code_map[res]['codons'][se[i] % self._get_nc(res)] for i, res in enumerate(self.config.seq)]) #HERE
        n_seqs = [get_seq(se) for se in members]
        return n_seqs
