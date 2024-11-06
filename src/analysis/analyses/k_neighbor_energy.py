<<<<<<< HEAD
import copy
import itertools
from src.analysis.analysis import Analysis
from src.rna_folding.rna_folder import RNAFolder
from src.rna_structure.structure import RNAStructure
from src.rna_structure.structure_io import StructureIO
from src.rna_structure.structure_convert import StructureConvert
=======
from src.analysis.analysis import Analysis
>>>>>>> bfe396c (added analysis)


class kNeighborEnergySearch(Analysis):
    """
    For a given input sequence and structure (connectivity table or TODO dot-bracket),
    calculate distribution of energies by changing k neighbors (stems)

    Parameters
    ----------
    config : AnalysisParser
        Object containing user inputs

    """

    def __init__(self, config):
        super().__init__(config)
        # read sequence from structure file
        self.connect_table = StructureIO()._ct_to_dataframe(self.config.args.input)
        seq = "".join([self.connect_table["Nucleotide"].iloc[i] for i in range(self.connect_table["Index"].iloc[-1])])

        # generate stems for sequence - need min stem length and min loop length!
        self.rna_folder_obj = RNAFolder(config)
        self.rna_folder_obj._declare_stem_vars(seq)

        # convert from connectivity table to stem tuples
        self.rna_struct_obj = RNAStructure()
        self.struct_conv_obj = StructureConvert()
        stems = self.struct_conv_obj._connect_table_to_stems(
            self.rna_folder_obj.n, self.connect_table
        )

        self.rna_folder_obj._gen_stems()

        # find indices of stems from connect table and place in sorted list
        self.active_stem_indices = []
        for stem in stems:
            found = False
            for j in range(len(self.rna_folder_obj.stems)):
                if stem == self.rna_folder_obj.stems[j]:
                    found = True
                    self.active_stem_indices.append(j)
            # if stem in connect table is not found in possible stem list
            # (self.rna_folder_obj.stems), likely because value for -ms is
            # higher than some stem lengths in input structure or a loop
            # smaller than specified by -ml is observed. Then add
            # observed stem to stem list and active stem list. This way,
            # the observed stems will be considered in the bit flipping
            # process
            if not found:
                self.rna_folder_obj.stems.append(stem)
                self.active_stem_indices.append(len(self.rna_folder_obj.stems)-1)
        self.active_stem_indices.sort()
        print(self.rna_folder_obj.stems)
        print(self.active_stem_indices)

        self.rna_folder_obj.len_stem_list = len(self.rna_folder_obj.stems)
        self.rna_folder_obj._compute_h_and_J()

        # keep list of calculated energies, starting with the initial structure
        self.energies = [self.rna_folder_obj._calc_score(self.active_stem_indices)]

        self._analyze()

    def _analyze(self):
        # brute force solution for k = 1
        new_stems = []
        valid_neighbor_count = 0
        lower_energy_neighbor_count = 0
        for i in range(self.rna_folder_obj.len_stem_list):
            if i in self.active_stem_indices:
                new_stems = [j for j in self.active_stem_indices if j != i]
            elif i not in self.active_stem_indices:
                new_stems = self.active_stem_indices + [i]

            new_energy = self.rna_folder_obj._calc_score(new_stems)
            if new_energy < 0:
                valid_neighbor_count += 1
            if new_energy < self.energies[0]:
                lower_energy_neighbor_count += 1
            self.energies.append(new_energy)
        #print(self.energies)
        print(self.energies[0], valid_neighbor_count, lower_energy_neighbor_count, min(self.energies))
