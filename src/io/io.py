class IO(object):
    """
    Class containing logic for reading and writing specific output file
    types

    """

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
