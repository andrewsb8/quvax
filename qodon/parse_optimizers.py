from qodon.optimizers.tf_differential_evo import TfDiffEv

def parse_optimizers(object):
    if object.args.codon_optimizer == "TFDE":
        TfDiffEv(object)
    else:
        print("Error here.")
