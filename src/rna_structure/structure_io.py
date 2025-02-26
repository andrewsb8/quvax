import pandas as pd


class StructureIO(object):
    """
    Class containing logic for reading and writing specific output file
    types

    """

    def _read_dot_bracket(self, dot_bracket_file):
        """
        Reads a dot bracket file which has lines
        - comment
        - sequence
        - structure in dot bracket notation

        Returns the second and third lines

        """
        dot_bracket = open(dot_bracket_file, 'r')
        data = dot_bracket.readlines()
        dot_bracket.close()
        seq = data[-2]
        struct = data[-1]
        return seq, struct

    def _ct_to_dataframe(self, ct_file):
        """
        Converts connectivity table from file to dataframe

        """
        df = pd.read_csv(ct_file, delim_whitespace=True, skiprows=1, header=None)
        df.columns = [
            "Index",
            "Nucleotide",
            "Previous",
            "Next",
            "Paired With",
            "Counter",
        ]
        return df

    def _get_sequence_from_connect_table(self, ct_dataframe):
        """
        Returns sequence from a dataframe collected from a connectivity
        table. Can only be run after _ct_to_dataframe!

        """
        return "".join(
            [
                ct_dataframe["Nucleotide"].iloc[i]
                for i in range(ct_dataframe["Index"].iloc[-1])
            ]
        )

    def _write_dot_bracket(self, output_file, energy, sequence, structure):
        """
        Write a simple dot bracket file containing the energy, sequence,
        and secondary structure in dot bracket notation

        """
        output = open(output_file, "w")
        output.write("> Folded energy: " + str(energy) + "\n")
        output.write(sequence + "\n")
        output.write(structure + "\n")
        output.close()

    def _write_connect_table(self, output_file, sequence, energy, pairs):
        """
        Write a connectivity table file.

        pairs : list
            list where the index+1 refers to base sequence number and the
            value at that index indicates the base pair sequence number.
            zero refers to unpaired base.

        """
        output = open(output_file, "w")
        output.write(f"{len(sequence):5} Energy: {energy} output from QuVax\n")
        for i in range(len(pairs)):
            output.write(f"{i+1:5} {sequence[i]} {i:5} {i+2:5} {pairs[i]:5} {i+1:5}\n")
        output.close()
