from src.analysis.analysis import Analysis
from src.rna_folding.rna_folder import RNAFolder


class StemSaturation(Analysis):
    """
    Reads an RNA structure (i.e. connectivity table) and enumerates
    the number of stems which can generated from the set of unpaired
    bases. Paired bases are replaced in the sequence string by an 'X'
    and then all possible stems are enumerated given the modified sequence,
    where 'X's will not be considered. Structures which cannot have any
    more stems added are considered "saturated".

    Parameters
    ----------
    config : AnalysisParser
        Object containing user inputs

    """

    def __init__(self, config):
        super().__init__(config)
        self.rna_folder = RNAFolder(config)
        structure = self.rna_folder._ct_to_dataframe(self.config.args.input)
        self.base_pairs = self.rna_folder._get_base_pair_indices(structure)
        self.seq = self.rna_folder._get_sequence_from_connect_table(structure)
        self._analyze()

    def _analyze(self):
        self.config.log.info("Input sequence: " + self.seq)
        list_seq = list(self.seq) # need to do this because strings are immutable
        for i in range(len(self.seq)):
            if self.base_pairs[i] != 0:
                list_seq[i] = 'X' # this character will be ignored in stem generation
        self.seq = "".join(list_seq)
        self.rna_folder._declare_stem_vars(self.seq)
        self.rna_folder._gen_stems()
        self.n_stems = len(self.rna_folder.stems)
        self.config.log.info("Modified sequence: " + self.seq)
        self.config.log.debug("Possible stems to be added: " + str(self.rna_folder.stems))
        self.config.log.info("Number of stems can be added to structure:")
        self.config.log.debug(str(self.n_stems))
        print(self.n_stems, flush=True)
