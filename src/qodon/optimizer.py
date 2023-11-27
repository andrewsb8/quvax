from abc import ABC, abstractmethod
from src.params.parser import Parser
from src.rna_folding.rna_folders.simulated_annealer import QuantumSimAnnealer
import python_codon_tables as pct
from Bio.Seq import Seq
import random
import numpy as np
import pandas as pd
from typing import List
import pickle


class CodonOptimizer(ABC):
    """
    Parent class for all codon optimizer classes.

    Parameters
    ----------
    config : Parser
        Object containing user inputs
    code_map : List
        Map of amino acids to codons based on species
    initial_sequences : List
        Randomly generated codon sequences for input amino acid sequence based on code map

    """
    def __init__(self, config: Parser):
        self.config = config
        self.config.log.info("Beginning codon optimization")
        self.codon_table, self.codon_scores, self.code_map = self._construct_codon_table()
        self.initial_sequences = self._generate_sequences(self.config.args.n_trials)
        self.optimization_process = {'generation-size': self.config.args.n_trials,
                                     'optimizer':       self.config.args.codon_optimizer,
                                     'sequences':       [], #nested list of sequences
                                     'scores':          [], #list of scores where index corresponds to sequence
                                     'sec_struct':      []} #list of secondary structure information for a sequence

    @abstractmethod
    def _optimize(self):
        pass

    def _convert_to_p_list(self, a):
        '''
        Helper function

        '''
        j = 0
        l = []
        for i in a:
            j += i
            l.append(j)
        return l


    def _construct_codon_table(self):
        '''
        Build reference table containing:

            amino acid:codon mappings
            codon frequencies

        This data is referenced by both the GA and the BQM

        '''
        # Load codon data
        codons_tables = pct.get_all_available_codons_tables()
        table = pct.get_codons_table(self.config.args.species)
        df = pd.DataFrame([(a, c, s) for a, v in table.items()
                           for c, s in v.items() if a != '*'],
                          columns=['aa', 'codon', 'score'])

        # Transform data into useful format
        df['tup'] = df.apply(lambda x: (x['codon'], x['score']), axis=1)
        by_aa = df.groupby('aa')
        ms_by_aa = by_aa['tup'].apply(list).apply(
            lambda x: max(x, key=lambda l: l[1]))
        df['log_score'] = df.apply(
            lambda x: abs(np.log(x['score'] / ms_by_aa.loc[x['aa']][1])), axis=1)

        # Merge lists of data into dataframe
        code_map_2 = pd.DataFrame(by_aa['score'].apply(list))
        code_map_2 = code_map_2.merge(pd.DataFrame(by_aa['codon'].apply(list)),
                                      left_index=True,
                                      right_index=True)
        code_map_2 = code_map_2.merge(pd.DataFrame(by_aa['log_score'].apply(list)),
                                      left_index=True,
                                      right_index=True)
        code_map_2.rename(columns={
            'score': 'scores',
            'codon': 'codons',
            'log_score': 'log_scores'
        },
                          inplace=True)

        # Convert dataframe to dict for quick lookups
        code_map_2['probs'] = code_map_2['scores'].apply(self._convert_to_p_list)
        code_map = code_map_2.to_dict('index')
        codon_scores = dict([
            item for sublist in
            [list(zip(_['codons'], _['log_scores'])) for _ in code_map.values()]
            for item in sublist
        ])
        codon_scores = {k: abs(v) for k, v in codon_scores.items()}

        codon_table = {k: v['codons'] for k, v in code_map.items()}

        return codon_table, codon_scores, code_map

    def _generate_sequences(self, ntrials) -> List:
        initial_members = []
        for i in range(ntrials):
            d_sequence = ""
            chosen_indices = []
            for res in self.config.seq:
                random_prob = random.uniform(0.0, 1.0)
                reference_chances = self.code_map[res]['probs']
                passing_indices = []
                for chance in reference_chances:
                    if chance > random_prob:
                        passing_indices.append(reference_chances.index(chance))
                chosen_index = passing_indices[0]
                chosen_indices.append(chosen_index)
                d_sequence += self.code_map[res]['codons'][chosen_index]
            member = chosen_indices
            initial_members.append(member)
        return initial_members

    def _fold_rna(self, nseq):
        '''
        Compute Minimum Free Energy (MFE) of RNA fold.

        '''
        folded_rna = QuantumSimAnnealer(nseq, self.config)
        return folded_rna.best_score

    def _extend_output(self, sequences, scores, sec_struct):
        self.optimization_process['sequences'].extend(sequences)
        self.optimization_process['scores'].extend(scores)
        return

    def _pickle_output(self):
        output_file = open(self.config.args.output, "wb")
        pickle.dump(self.optimization_process, output_file)
        return

    def _read_pickle(self):
        #read previous optimization and continue process. not ready yet
        return

    def _get_num_codons(self, res):
        '''
        Extract number of possible codons for each amino acid

        '''
        return len(self.code_map[res]['codons'])

    def _reverse_translate(self, sequence):
        '''
        Convert to nucleotide sequence from integer indices of code map

        '''

        return ''.join([self.code_map[res]['codons'][sequence[i] % self._get_num_codons(res)] for i, res in enumerate(self.config.seq)])

    def _verify_dna(self, sequence):
        '''
        Translate nucleotide sequence to make sure it matches input

        '''
        if self.config.seq != str(Seq(sequence).transcribe().translate()):
            self.config.log.error("Error: Codon sequence did not translate properly!")
            raise ValueError(
                "Error: Codon sequence did not translate properly!")
        else:
            self.config.log.info("Final codon sequence translated properly.")
            return True

    def _get_optimized_sequence(self):
        '''
        get lowest energy and associated sequence from all sequences generated

        '''
        self.mfe = np.min(self.optimization_process['scores'])
        self.mfe_index = np.argmin(self.optimization_process['scores'])
        self.final_codons = self.optimization_process['sequences'][self.mfe_index]
        if self._verify_dna(self.final_codons):
            self.config.log.info("Minimum energy codon sequence: " + self.final_codons)
            self.config.log.info("Energy of codon sequence: " + str(self.mfe))
