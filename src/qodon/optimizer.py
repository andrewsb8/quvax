from src.params.design_parser import DesignParser
from src.rna_folding.rna_folders.simulated_annealer import SimulatedAnnealer
from src.rna_folding.rna_folders.classical_mc import MC
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
        self.folder = self._determine_folder()
        if not self.config.args.resume:
            if self.config.args.target is not None:
                self._verify_target()
                self._fold_target()
            self.initial_sequences = self._generate_sequences(self.config.args.n_trials)
            self.mfe = 1000000  # set min free energy to high number
        else:
            self._load_random_state()
            self.mfe = self.config.mfe
            self.initial_sequences = self.config.initial_sequences
            if self.config.args.target is not None:
                self.target_folded_energy = self.config.target_folded_energy

    @abstractmethod
    def _optimize(self):
        pass

    def _update_codon_step(self):
        self.codon_optimize_step += 1
        if self.config.args.resume:
            step = self.codon_optimize_step + self.config.generations_sampled
            total = self.config.args.codon_iterations + self.config.generations_sampled
        else:
            step = self.codon_optimize_step
            total = self.config.args.codon_iterations
        sys.stderr.write(
            "Codon optimization step ("
            + str(step)
            + ") of total steps ("
            + str(total)
            + ")\r"
        )
        if self.codon_optimize_step == self.config.args.codon_iterations:
            sys.stderr.write("\n")
            self.config.log.info(
                "Number of Generations: ("
                + str(step)
                + ") of total generations ("
                + str(total)
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
            chosen_indices = []
            for res in self.config.seq:
                chosen_indices.append(
                    random.randint(0.0, len(self.code_map[res]["codons"]) - 1)
                )
            initial_members.append(chosen_indices)
        return initial_members

    def _determine_folder(self):
        if self.config.args.solver == "SA":
            return SimulatedAnnealer(self.config)
        elif self.config.args.solver == "MC":
            return MC(self.config)
        else:
            raise NotImplementedError("Invalid RNA folding solver. Use -h for available options.")

    def _fold_rna(self, nseq):
        """
        Compute Minimum Free Energy (MFE) of RNA fold.

        """
        self.folder._fold(nseq)
        return self.folder.best_score

    def _write_output(self, sequences, energies, secondary_structure):
        if self.config.args.resume:
            step = self.codon_optimize_step + self.config.generations_sampled
        else:
            step = self.codon_optimize_step
        for i in range(len(energies)):
            self.config.db_cursor.execute(
                "INSERT INTO OUTPUTS(sim_key, population_key, generation, sequences, energies) VALUES(?, ?, ?, ?, ?);",
                (
                    self.config.sim_key,
                    i,
                    step,
                    sequences[i],
                    energies[i],
                ),
            )
            self.config.db.commit()
        return

    def _get_num_codons(self, res):
        """
        Extract number of possible codons for each amino acid

        """
        return len(self.code_map[res]["codons"])

    def _convert_ints_to_codons(self, sequence):
        """
        Convert to nucleotide sequence from integer indices of code map

        """

        return "".join(
            [
                self.code_map[res]["codons"][sequence[i] % self._get_num_codons(res)]
                for i, res in enumerate(self.config.seq)
            ]
        )

    def _convert_codons_to_ints(self, sequence):
        """
        Convert to integer indices from nucleotide sequence of code map

        """

        return [
            self.code_map[res]["codons"].index(sequence[i * 3 : (i * 3) + 3])
            for i, res in enumerate(self.config.seq)
        ]

    def _iterate(self, sequences):
        """
        Function containing references to the steps taken in each codon
        optimization iteration: convert codon integer sequences to codon
        strings, calculate the folding energies of the sequences, writes each
        to the database, and updates the min free energy.

        """
        self.list_seqs = [self._convert_ints_to_codons(s) for s in sequences]
        self.energies = [self._fold_rna(s) for s in self.list_seqs]
        self._update_mfe(self.energies)
        self._write_output(self.list_seqs, self.energies, None)

    def _post_process(self):
        if self.config.args.resume:
            num = self.codon_optimize_step + self.config.generations_sampled
            # below is the equivalent of TRUNCATE in MySQL DB, SQLite has
            # different syntax. Clear table to avoid duplicate degenerate
            # sequences
            self.config.db_cursor.execute("DELETE FROM MFE_SEQUENCES;")
        else:
            # add one to account for initial sequences
            num = self.codon_optimize_step
        self.config.db_cursor.execute(
            "UPDATE SIM_DETAILS SET generations_sampled = ? WHERE protein_sequence = ?;",
            (num, self.config.seq),
        )
        self.config.db.commit()
        self._save_random_state()
        self._get_number_unique_sequences()
        self._get_optimized_sequences()
        if self.config.args.target is not None:
            self._check_target()

    def _load_random_state(self):
        """
        Funciton to restore the random number generator to the state it was in
        at the end of the previous execution of design.py. The file name was
        retrieved from the database containing the simulation details.

        """
        file = open(self.config.args.state_file, "rb")
        state = pickle.load(file)
        random.setstate(state)
        file.close()

    def _save_random_state(self):
        file = open(self.config.args.state_file, "wb")
        pickle.dump(random.getstate(), file)
        file.close()

    def _get_number_unique_sequences(self):
        if self.config.args.resume:
            step = self.codon_optimize_step + self.config.generations_sampled
        else:
            step = self.codon_optimize_step
        self.config.db_cursor.execute("SELECT COUNT(DISTINCT sequences) from OUTPUTS;")
        num = self.config.db_cursor.fetchall()[0][0]
        # step+1 to account for initial randomly generated sequences
        self.config.log.info(
            "Number of unique sequences sampled: "
            + str(num)
            + " of possible "
            + str((step + 1) * self.config.args.n_trials)
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
