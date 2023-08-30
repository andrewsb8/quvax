from qodon.optimizers.tf_differential_evo import TfDiffEv

def parse_optimizers(optimizer):
    if optimizer == "TFDE":
        TfDiffEv()
    else:
        print("Error here.")
