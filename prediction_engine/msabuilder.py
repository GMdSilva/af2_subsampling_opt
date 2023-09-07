"""
Defines class for building MSA for subsampled AF2 predictions using jackhmmer.
"""

import os

from user_settings.config import JACKHMMER_PATH


class MSABuilder:
    """
    Class for building MSA for subsampled AF2 predictions using jackhmmer.
    """
    def __init__(self, prefix: str, sequence: str = 'HLDERPLACE'):
        """
        Args:
        - prefix (str): Usually the name of the target protein.
        - sequence (list): Target sequence.
        """
        self.prefix = prefix
        self.sequence = sequence

    def _reformat_sequences(self, msa: list) -> list:
        """
        Format sequences in the MSA.

        Args:
            msa (list): Multi-sequence alignment data.

        Returns:
            list: Formatted sequences.
        """
        formatted_sequences = []
        for idx, sequence in enumerate(msa[0]):
            formatted_sequence = f">sequence_{idx + 1}\n{sequence}\n"
            formatted_sequences.append(formatted_sequence)
        return formatted_sequences

    def jackhmmer_to_fasta(self, msa: list):
        """
        Convert the MSA data from jackhmmer to fasta format.

        Args:
            msa (list): Multi-sequence alignment data.
        """
        results_dir = os.path.join("results",
                                   "msas",
                                   f"{self.prefix}_hmmer_{len(msa[0])}_seqs.a3m")
        formatted_sequences = self._reformat_sequences(msa)

        with open(results_dir, 'w', encoding="utf-8") as file:
            for sequence in formatted_sequences:
                file.write(sequence)

    def run_jackhmmer(self):
        """
        Execute the jackhmmer command. (Implementation is not provided)
        """
        # NOT IMPLEMENTED #
        # msa = run_hmmer(self.sequence)
        msa = [['PLACEHOLDER', 'HOLDERPLACE']]
        self.jackhmmer_to_fasta(msa)

    def build_jackhmmer_msa(self):
        """
        Build an MSA for a target sequence using jackhmmer
            (querying UniRef90, small BFD, and mgnify remotely).
        """
        jackhmmer_installed = os.path.isfile(JACKHMMER_PATH)

        if not jackhmmer_installed:
            print('JACKHMMER not found, not building MSA')
        else:
            self.run_jackhmmer()
