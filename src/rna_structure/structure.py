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
