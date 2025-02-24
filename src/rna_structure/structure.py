import numpy as np


class RNAStructure(object):
    """
    Class containing methods regarding RNA secondary structure

    """

    def _get_all_interactions(self):
        return self._get_wc_interactions() + self._get_wobble_interactions()

    def _get_wc_interactions(self):
        return [
            ("A", "U"),
            ("U", "A"),
            ("G", "C"),
            ("C", "G"),
        ]

    def _get_wobble_interactions(self):
        return [
            ("G", "U"),
            ("U", "G"),
        ]

    def _calc_stem_separation(self, first_base, last_base, stem_length):
        """
        Calculates the number of bases between the closest bases in a stem.
        The below uses stem_length - 1 because of the following:

        Say we have a stem tuple (1, 8, 3) representing the base pairs
        (1,8), (2,7), (3,6). The stem is of length 3, but 3-1 is 2.

        So, in (i, j, k) notation, to compare the distance between base
        3 and six I need to do (8 - (3-1)) - (1 + (3-1)) = 6-3 = 3. Then,
        need to subtract again by 1 because only bases 4 and 5 are between
        all of the bases in the stem, so the distance is then 3-1 = 2 bases.

        The value can be used for check loop sizes in hairpins and also for
        enforcing span requirements.

        """
        return (last_base - (stem_length - 1)) - (first_base + (stem_length - 1)) - 1

    def _get_base_pairs(self, structure):
        """
        Return list of sequence location of base pair from
        pandas dataframe ```structure```. Index of list is position
        in sequence - 1.

        """
        return structure["Paired With"]

    def _detect_stem_overlap(self, stem1, stem2):
        """
        Function to detect if a base is involved in two stems which is
        a stem overlap. Returns True if the stems overlap and False
        otherwise.

        """
        from src.rna_structure.structure_convert import StructureConvert

        struct_convert = StructureConvert()
        pairs1 = struct_convert._stem_to_pair_list(stem1)
        pairs2 = struct_convert._stem_to_pair_list(stem2)

        for pair1 in pairs1:
            for pair2 in pairs2:
                if pair1[0] in pair2 or pair1[1] in pair2:
                    return True
        return False

    def _is_pseudo(self, stem1, stem2):
        first = np.argmin([stem1[0], stem2[0]])
        second = np.argmax([stem1[0], stem2[0]])
        stem_pair = [stem1, stem2]
        if (
            stem_pair[first][0] <= stem_pair[second][0] <= stem_pair[first][1]
            and stem_pair[second][0] <= stem_pair[first][1] <= stem_pair[second][1]
        ):
            return True
        return False
