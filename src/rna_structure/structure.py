import numpy as np


class RNAStructure(object):
    """
    Class containing methods regarding RNA secondary structure

    """

    def _get_all_interactions(self):
        print(self._get_wc_interactions() + self._get_wobble_interactions())
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

    def _detect_stem_of_loops(self, ref_stem, stems):
        """
        Function to determin if stem (pseudoknot) if formed
        by bases which are both in hairpin loops

        """

        stem_in_loop_count = 0
        for stem in stems:
            if stem == ref_stem:
                continue
            # check if number of bases between base pairs is
            # below arbitrary threshold which indicates a tight,
            # inflexible loop which should not be able to form
            # sufficient stem lengths with other loops
            if (stem[1] - stem[2]) - (stem[0] + stem[2]) <= 15:
                # then check to see if ref stem falls between it
                if (ref_stem[0] > stem[0] + (stem[2]-1) and ref_stem[0] + (ref_stem[2]-1) < stem[1] - (stem[2]-1)): # or (ref_stem[1] - (ref_stem[2] -1 ) > stem[0] + (stem[2]-1) and ref_stem[1] < stem[1] - (stem[2]-1)):
                    stem_in_loop_count += 1
                elif (ref_stem[1] - (ref_stem[2] - 1) > stem[0] + (stem[2]-1) and ref_stem[1] < stem[1] - (stem[2] -1)):
                    stem_in_loop_count += 1
            if stem_in_loop_count > 1:
                    return True
        return False
