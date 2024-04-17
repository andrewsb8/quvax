from src.params.design_parser import DesignParser
from src.rna_folding.rna_folders.simulated_annealer import SimulatedAnnealer
from abc import ABC, abstractmethod
import python_codon_tables as pct
from Bio.Seq import Seq
import random
import numpy as np
import pandas as pd
from typing import List
import pickle
import sys


class CodonOptimizer(ABC):
    """
    Parent class for all codon optimizer classes.

    Parameters
    ----------
    config : DesignParser
        Object containing user inputs
    code_map : List
        Map of amino acids to codons based on species
    initial_sequences : List
        Randomly generated codon sequences for input amino acid sequence based on code map

    """

    def __init__(self, config: DesignParser):
        self.config = config
        random.seed(self.config.args.random_seed)
        self.codon_optimize_step = 0
        self.config.log.info("Beginning codon optimization")
        (
            self.codon_table,
            self.codon_scores,
            self.code_map,
        ) = self._construct_codon_table()
        if self.config.args.target is not None:
            self._verify_target()
            self._fold_target()
        if self.config.args.resume == False:
            self.initial_sequences = self._generate_sequences(self.config.args.n_trials)
        else:
            self.initial_sequences = self.config.initial_sequences
        self.mfe = 1000000  # set min free energy to high number

    @abstractmethod
    def _optimize(self):
        pass

    def _update_codon_step(self):
        self.codon_optimize_step += 1
        sys.stderr.write(
            "Codon optimization step ("
            + str(self.codon_optimize_step)
            + ") of total steps ("
            + str(self.config.args.codon_iterations)
            + ")\r"
        )
        if self.codon_optimize_step == self.config.args.codon_iterations:
            sys.stderr.write("\n")
            self.config.log.info(
                "Number of Generations: ("
                + str(self.codon_optimize_step)
                + ") of total generations ("
                + str(self.config.args.codon_iterations)
                + ")\n"
            )

    def _update_mfe(self, energies):
        for energy in energies:
            if energy < self.mfe:
                self.mfe = energy

    def _convert_to_p_list(self, a):
        """
        Helper function

        """
        j = 0
        l = []
        for i in a:
            j += i
            l.append(j)
        return l

    def _construct_codon_table(self):
        """
        Build reference table containing amino acid->codon mappings, codon frequencies

        """
        # Load codon data
        # codons_tables = pct.get_all_available_codons_tables()
        table = pct.get_codons_table(self.config.args.species)
        df = pd.DataFrame(
            [(a, c, s) for a, v in table.items() for c, s in v.items() if a != "*"],
            columns=["aa", "codon", "score"],
        )

        # Transform data into useful format
        df["tup"] = df.apply(lambda x: (x["codon"], x["score"]), axis=1)
        by_aa = df.groupby("aa")
        ms_by_aa = by_aa["tup"].apply(list).apply(lambda x: max(x, key=lambda l: l[1]))
        df["log_score"] = df.apply(
            lambda x: abs(np.log(x["score"] / ms_by_aa.loc[x["aa"]][1])), axis=1
        )

        # Merge lists of data into dataframe
        code_map_2 = pd.DataFrame(by_aa["score"].apply(list))
        code_map_2 = code_map_2.merge(
            pd.DataFrame(by_aa["codon"].apply(list)), left_index=True, right_index=True
        )
        code_map_2 = code_map_2.merge(
            pd.DataFrame(by_aa["log_score"].apply(list)),
            left_index=True,
            right_index=True,
        )
        code_map_2.rename(
            columns={"score": "scores", "codon": "codons", "log_score": "log_scores"},
            inplace=True,
        )

        # Convert dataframe to dict for quick lookups
        code_map_2["probs"] = code_map_2["scores"].apply(self._convert_to_p_list)
        code_map = code_map_2.to_dict("index")
        codon_scores = dict(
            [
                item
                for sublist in [
                    list(zip(_["codons"], _["log_scores"])) for _ in code_map.values()
                ]
                for item in sublist
            ]
        )
        codon_scores = {k: abs(v) for k, v in codon_scores.items()}

        codon_table = {k: v["codons"] for k, v in code_map.items()}

        return codon_table, codon_scores, code_map

    def _generate_sequences(self, ntrials) -> List:
        initial_members = []
        for i in range(ntrials):
            d_sequence = ""
            chosen_indices = []
            for res in self.config.seq:
                random_prob = random.uniform(0.0, 1.0)
                reference_chances = self.code_map[res]["probs"]
                passing_indices = []
                for chance in reference_chances:
                    if chance > random_prob:
                        passing_indices.append(reference_chances.index(chance))
                chosen_index = passing_indices[0]
                chosen_indices.append(chosen_index)
                d_sequence += self.code_map[res]["codons"][chosen_index]
            member = chosen_indices
            initial_members.append(member)
        return initial_members

    def _fold_rna(self, nseq):
        """
        Compute Minimum Free Energy (MFE) of RNA fold.

        """
        folded_rna = SimulatedAnnealer(nseq, self.config)
        return folded_rna.best_score

    def _write_output(self, sequences, energies, secondary_structure):
        for i in range(len(energies)):
            self.config.db_cursor.execute(
                "INSERT INTO OUTPUTS(sim_key, population_key, generation, sequences, energies) VALUES(?, ?, ?, ?, ?);",
                (
                    self.config.sim_key,
                    i,
                    self.codon_optimize_step,
                    sequences[i],
                    energies[i],
                ),
            )
            self.config.db.commit()
        return

    def _read_output(self):
        # read previous optimization and continue process.
        raise NotImplementedError()

    def _get_num_codons(self, res):
        """
        Extract number of possible codons for each amino acid

        """
        return len(self.code_map[res]["codons"])

    def _reverse_translate(self, sequence):
        """
        Convert to nucleotide sequence from integer indices of code map

        """

        return "".join(
            [
                self.code_map[res]["codons"][sequence[i] % self._get_num_codons(res)]
                for i, res in enumerate(self.config.seq)
            ]
        )

    def _verify_dna(self, sequence):
        """
        Translate nucleotide sequence to make sure it matches input

        """
        if self.config.seq != str(Seq(sequence).transcribe().translate()):
            self.config.log.error(
                "Error: Codon sequence did not translate properly! Sequence: "
                + sequence
            )
        else:
            self.config.log.info("Codon sequence translated properly.")

    def _get_optimized_sequences(self):
        """
        Get lowest energy sequences from all sampled sequences

        """
        # write min free energy to log and db
        self.config.log.info("Minimum energy of codon sequences: " + str(self.mfe))
        self.config.db_cursor.execute(
            "UPDATE SIM_DETAILS SET min_free_energy = ? WHERE protein_sequence = ?;",
            (self.mfe, self.config.seq),
        )
        self.config.db.commit()

        # get number and list of degenerate min free energy sequences
        self.config.db_cursor.execute(
            f"SELECT COUNT(sequences) FROM OUTPUTS WHERE energies = {self.mfe};"
        )
        num_degen_sequences = self.config.db_cursor.fetchall()[0][0]
        self.config.log.info(
            "Number of degenerate minimum free energy sequences sampled: "
            + str(num_degen_sequences)
        )
        self.config.db_cursor.execute(
            "INSERT INTO MFE_SEQUENCES (sequences) SELECT sequences FROM OUTPUTS WHERE energies = ?",
            (self.mfe,),
        )
        self.config.db.commit()
        self.config.log.info("Finished parsing optimized sequences.")

    def _verify_target(self):
        """
        Verify target codes for the input protein sequence

        """
        self.config.log.info("Verifying target codes for input protein.")
        self._verify_dna(self.config.args.target)

    def _fold_target(self):
        self.target_folded_energy = self._fold_rna(self.config.args.target)
        self.config.log.info(
            "Target sequence folding energy: " + str(self.target_folded_energy)
        )
        self.config.db_cursor.execute(
            "UPDATE SIM_DETAILS SET target_min_free_energy = ? WHERE target_sequence = ?;",
            (self.target_folded_energy, self.config.args.target),
        )
        self.config.db.commit()
        self.config.log.info("\n")

    def _check_target(self):
        """
        Check if target codon sequence was sampled and if it was lowest energy

        """
        self.config.db_cursor.execute(
            f"SELECT COUNT(sequences) FROM MFE_SEQUENCES WHERE sequences = '{self.config.args.target}';"
        )
        mfe_samples = self.config.db_cursor.fetchall()[0][0]
        if mfe_samples > 0:
            self.config.log.info(
                "The target codon sequence is in the list of minimum free energy sequences!"
            )
            return  # return early if condition is met to avoid unnecessary query

        self.config.db_cursor.execute(
            f"SELECT COUNT(sequences) FROM OUTPUTS WHERE sequences = '{self.config.args.target}';"
        )
        samples = self.config.db_cursor.fetchall()[0][0]
        if mfe_samples == 0 and samples > 0:
            self.config.log.warning(
                "The target codon sequence was sampled but was not the lowest free energy sequence."
            )
        else:
            self.config.log.error("The target codon sequence was NOT sampled.")
