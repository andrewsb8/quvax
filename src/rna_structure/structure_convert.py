class StructureConvert(object):
    """
    Class containing logic for converting between different
    representations of RNA structure

    """

    def _stems_to_dot_bracket(self, sequence_len, stems):
        """
        Function to convert a list of stems in a sequence to a dot-bracket
        notation.

        """

        dot_bracket = ["." for i in range(sequence_len)]
        for stem in stems:
            stem_pair_list = self._stem_to_pair_list(stem)
            for i in range(len(stem_pair_list)):
                dot_bracket[stem_pair_list[i][0] - 1] = "("
                dot_bracket[stem_pair_list[i][1] - 1] = ")"

        # check for pseudoknots. pseudoknot cannot be detected if there is only one stem
        for i in range(len(stems)):
            for j in range(i + 1, len(stems)):
                if self._is_pseudo(stems[i], stems[j]):
                    stem_pair_list = self._stem_to_pair_list(stems[j])
                    for k in range(len(stem_pair_list)):
                        dot_bracket[stem_pair_list[k][0] - 1] = "["
                        dot_bracket[stem_pair_list[k][1] - 1] = "]"

        return "".join(dot_bracket)

    def _stems_to_connect_list(self, sequence_len, stems):
        """
        Helper function to create a "connect list" from a list of
        stems.

        pair : list
            list where the index+1 refers to base sequence number and the
            value at that index indicates the base pair sequence number.
            zero refers to unpaired base. This is passed to class IO in
            src.io.io to help write connectivity table files.

        """
        connect_list = [0 for i in range(sequence_len)]
        for stem in stems:
            stem_pair_list = self._stem_to_pair_list(stem)
            for i in range(len(stem_pair_list)):
                connect_list[stem_pair_list[i][0] - 1] = stem_pair_list[i][1]
                connect_list[stem_pair_list[i][1] - 1] = stem_pair_list[i][0]
        return connect_list
