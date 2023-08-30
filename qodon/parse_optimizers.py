from design import QuDesign
from qodon.optimizers.tf_differential_evo import TfDiffEv

class ParseOptimizers(QuDesign):
    def __init__(self):
        print(QuDesign.seq)
        if self.args.codon_optimizer == "TFDE":
            TfDiffEv()
        else:
            print("Error here.")
