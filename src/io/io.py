class IO(object):
    def _write_dot_bracket(self, output_file, energy, sequence, structure):
        output = open(output_file, "w")
        output.write("> Folded energy: " + str(energy) + "\n")
        output.write(sequence + "\n")
        output.write(structure + "\n")
        output.close()

    def _write_connect_table(self, output_file, sequence, pairs):
        # pairs is a list, see rna_folder for definition
        output = open(output_file, "w")
        output.write(f"{len(sequence):5} output from QuVax\n")
        for i in range(len(pairs)):
            output.write(f"{i+1:5} {sequence[i]} {i:5} {i+2:5} {pairs[i]:5} {i+1:5}\n")
        output.close()
