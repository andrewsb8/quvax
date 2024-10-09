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

    @staticmethod
    def _stem_to_pair_list(stem):
        """
        From stem definiton (see rna_folder.py), can generate list of tuples with
        indices of base pairs
        Ex:
        - self.stems[0] -> (1, 13, 3)
        - _stem_to_pair_list(self.stems[0]) -> [(1, 13), (2, 12), (3, 11)]
        - Above shows the list of indices of three base pairs comprising the stem

        """
        pair_list = []
        for ci in range(stem[2]):
            pair_list.append((stem[0] + ci, stem[1] - ci))
        return pair_list

    def _connect_table_to_stems(self, sequence_len, connect_table_df):
        stems = []
        i = 0
        while i < sequence_len:
            # only need to parse half of the dataframe so if detect a base pair index lower
            # than current index, can break the loop
            if connect_table_df["Paired With"].iloc[i] < connect_table_df["Index"].iloc[i]:
                break
            # if base pair is found, need to define the length of stem
            elif connect_table_df["Paired With"].iloc[i] != 0:
                # loop through connect table until nonsequential base pair is
                # found and then append to stems
                for j in range(i+1, sequence_len-1):
                    # if nonsequential base pair ordering or a unpaired base is found, use information to generate stem
                    if connect_table_df["Paired With"].iloc[j] != connect_table_df["Paired With"].iloc[j-1] - 1 or connect_table_df["Paired With"].iloc[j] == 0:
                        stems.append((connect_table_df["Index"].iloc[i], connect_table_df["Paired With"].iloc[i], j-i))
                        i += j
                        break
        return stems
