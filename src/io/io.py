class IO(object):
    def _write_dot_bracket(self, output_file, energy, sequence, structure):
        output = open(output_file, "w")
        output.write("> Folded energy: " + str(energy) + "\n")
        output.write(sequence + "\n")
        output.write(structure + "\n")
        output.close()
