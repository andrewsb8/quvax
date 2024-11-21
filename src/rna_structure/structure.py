import numpy as np


class RNAStructure(object):
    """
    Class containing methods regarding RNA secondary structure

    """

    def _get_wc_interactions(self):
        return [
            ("A", "U"),
            ("U", "A"),
            ("G", "C"),
            ("C", "G"),
            ("G", "U"),
            ("U", "G"),
        ]

    def _detect_stem_overlap(self, stem1, stem2):
        from src.rna_structure.structure_convert import StructureConvert

        struct_convert = StructureConvert()
        pairs1 = struct_convert._stem_to_pair_list(stem1)
        pairs2 = struct_convert._stem_to_pair_list(stem2)

        keep1 = [
            pair1
            for pair1 in pairs1
            if not any(_ in pair1 for _ in np.array(pairs2).flatten())
        ]
        keep2 = [
            pair2
            for pair2 in pairs2
            if not any(_ in pair2 for _ in np.array(pairs1).flatten())
        ]

        if len(keep1) == len(pairs1) and len(keep2) == len(pairs2):
            # No overlap
            return False
        else:
            # Overlap
            return True

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
