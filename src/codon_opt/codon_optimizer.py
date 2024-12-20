from src.config.config import Config
from abc import ABC, abstractmethod
from copy import deepcopy
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

    def __init__(self, config: Config):
        self.config = config
        random.seed(self.config.args.random_seed)
        self.config.log.debug(
            "Input protein sequence length: " + str(len(self.config.protein_sequence))
        )
        self.codon_optimize_step = 0
        self.convergence_count = 0
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
            self.initial_sequences = self._generate_sequences(
                self.config.args.population_size
            )
            self.members = deepcopy(self.initial_sequences)
            self.energies = [0 for member in self.members]
            self.sec_structs = ["" for member in self.members]
            self.min_free_energy = 1000000  # set min free energy to high number
            # tensorflow will calculate the energies for the initial population on its own
            if self.config.args.codon_optimizer != "TFDE":
                self._iterate(update_counter=False)
        else:
            self._load_random_state()
            self.convergence_count = self.config.args.convergence_count
            self.min_free_energy = self.config.min_free_energy
            self.initial_sequences = [
                self._convert_codons_to_ints(s) for s in self.config.initial_sequences
            ]
            self.members = deepcopy(self.initial_sequences)
            self.energies = self.config.energies
            self.sec_structs = self.config.sec_structs
            if self.config.args.target is not None:
                self.target_min_free_energy = self.config.target_min_free_energy
        self.config.log.debug("Beginning codon optimization")

    @abstractmethod
    def _optimize(self):
        pass

    def _update_codon_step(self, update_counter=True):
        if update_counter:
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
        new_min = False
        for energy in energies:
            if energy < self.min_free_energy:
                self.min_free_energy = energy
                # reset convergence_counter
                self.convergence_count = 0
                new_min = True
        if not new_min:
            self.convergence_count += 1

    def _check_convergence(self):
        if self.convergence_count >= self.config.args.convergence:
            self.config.log.info(
                "A new free minimum energy sequence has not been sampled in "
                + str(self.config.args.convergence)
                + " generations. Optimization is converged. Terminating.\n"
            )
            self._post_process()
            sys.exit(1)

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
        table = pct.get_codons_table(self.config.args.species, replace_U_by_T=False)
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
            for res in self.config.protein_sequence:
                chosen_indices.append(
                    random.randint(0.0, len(self.code_map[res]["codons"]) - 1)
                )
            initial_members.append(chosen_indices)
        return initial_members

    def _determine_folder(self):
        if self.config.args.solver == "SA":
            from src.rna_folding.rna_folders.simulated_annealer import SimulatedAnnealer

            return SimulatedAnnealer(self.config)
        elif self.config.args.solver == "MC":
            from src.rna_folding.rna_folders.classical_mc import MC

            return MC(self.config)
        elif self.config.args.solver == "ES":
            from src.rna_folding.rna_folders.exact_solver import ExactSolver

            return ExactSolver(self.config)
        else:
            raise NotImplementedError(
                "Invalid RNA folding solver. Use -h for available options."
            )

    def _fold_rna(self, nseq):
        """
        Compute Minimum Free Energy (MFE) of RNA fold.

        """
        self.folder._fold(nseq)

    def _write_output(self, sequences, energies, secondary_structure):
        if self.config.args.resume:
            step = self.codon_optimize_step + self.config.generations_sampled
        else:
            step = self.codon_optimize_step
        for i in range(len(energies)):
            self.config.db_cursor.execute(
                f"""INSERT INTO OUTPUTS(sim_key, population_key, generation,
                sequences, energies, secondary_structure) VALUES(
                '{self.config.sim_key}', '{i}', '{step}', '{sequences[i]}',
                '{energies[i]}', '{secondary_structure[i]}');"""
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
                for i, res in enumerate(self.config.protein_sequence)
            ]
        )

    def _convert_codons_to_ints(self, sequence):
        """
        Convert to integer indices from nucleotide sequence of code map

        """

        return [
            self.code_map[res]["codons"].index(sequence[i * 3 : (i * 3) + 3])
            for i, res in enumerate(self.config.protein_sequence)
        ]

    def _iterate(self, fold_sequences=True, update_counter=True):
        """
        Function containing references to the steps taken in each codon
        optimization iteration: convert codon integer sequences to codon
        strings, calculate the folding energies of the sequences, writes each
        to the database, and updates the min free energy.

        """

        self._update_codon_step(update_counter)
        self.list_seqs = [self._convert_ints_to_codons(s) for s in self.members]
        if fold_sequences:
            for s in range(len(self.list_seqs)):
                self.config.log.debug(
                    "Generation number: "
                    + str(self.codon_optimize_step)
                    + ". Population number "
                    + str(s)
                    + " of "
                    + str(self.config.args.population_size)
                )
                self._fold_rna(self.list_seqs[s])
                self.config.log.debug("Completed structure prediction.\n")
                self.energies[s] = self.folder.best_score
                self.sec_structs[s] = self.folder.dot_bracket
        self._update_mfe(self.energies)
        self._write_output(self.list_seqs, self.energies, self.sec_structs)
        if self.config.args.convergence > 0:
            self._check_convergence()
        if (
            self.codon_optimize_step != 0
            and self.config.args.checkpoint_interval != 0
            and self.codon_optimize_step % self.config.args.checkpoint_interval == 0
            and self.codon_optimize_step != self.config.args.codon_iterations
        ):
            sys.stderr.write("\n")
            self.config.log.info(
                "Writing checkpoint at step " + str(self.codon_optimize_step) + ":"
            )
            self._post_process()
            self.config.log.info("")

    def _post_process(self):
        if self.config.args.resume:
            num = self.codon_optimize_step + self.config.generations_sampled
            # below is the equivalent of TRUNCATE in MySQL DB, SQLite has
            # different syntax. Clear table to avoid duplicate degenerate
            # sequences
            self.config.db_cursor.execute(
                f"DELETE FROM MFE_SEQUENCES WHERE sim_key = '{self.config.sim_key}';"
            )
        else:
            # add one to account for initial sequences
            num = self.codon_optimize_step
        self.config.db_cursor.execute(
            f"UPDATE SIM_DETAILS SET generations_sampled = '{num}' WHERE sim_key = '{self.config.sim_key}';"
        )
        self.config.db_cursor.execute(
            f"UPDATE SIM_DETAILS SET convergence_count = '{self.convergence_count}' WHERE sim_key = '{self.config.sim_key}';"
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
        self.config.db_cursor.execute(
            f"SELECT COUNT(DISTINCT sequences) from OUTPUTS where sim_key = '{self.config.sim_key}';"
        )
        num = self.config.db_cursor.fetchall()[0][0]
        # "+ self.config.args.population_size" to account for initial randomly generated sequences
        self.config.log.info(
            "Number of unique sequences sampled: "
            + str(num)
            + " of possible "
            + str(
                (step * self.config.args.population_size)
                + self.config.args.population_size
            )
        )

    def _verify_dna(self, sequence):
        """
        Translate nucleotide sequence to make sure it matches input

        """
        if self.config.protein_sequence != str(Seq(sequence).transcribe().translate()):
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
        self.config.log.info(
            "Minimum energy of codon sequences: " + str(self.min_free_energy)
        )
        self.config.db_cursor.execute(
            f"UPDATE SIM_DETAILS SET min_free_energy = '{self.min_free_energy}' WHERE sim_key = '{self.config.sim_key}';"
        )
        self.config.db.commit()

        # get number and list of degenerate min free energy sequences
        self.config.db_cursor.execute(
            f"SELECT COUNT(sequences) FROM OUTPUTS WHERE energies = {self.min_free_energy} and sim_key = {self.config.sim_key};"
        )
        num_degen_sequences = self.config.db_cursor.fetchall()[0][0]
        self.config.log.info(
            "Number of degenerate minimum free energy sequences sampled: "
            + str(num_degen_sequences)
        )
        self.config.db_cursor.execute(
            f"INSERT INTO MFE_SEQUENCES (sim_key, sequences, secondary_structure) SELECT sim_key, sequences, secondary_structure FROM OUTPUTS WHERE energies = '{self.min_free_energy}'"
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
        self.target_min_free_energy = self._fold_rna(self.config.args.target)
        self.config.log.info(
            "Target sequence folding energy: " + str(self.target_min_free_energy)
        )
        self.config.db_cursor.execute(
            f"UPDATE SIM_DETAILS SET target_min_free_energy = '{self.target_min_free_energy}' WHERE target_sequence = '{self.config.args.target}';"
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
