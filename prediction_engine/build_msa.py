"""Functions for building MSA for subsampled AF2 predictions"""

import os
from user_settings.config import JACKHMMER_PATH, PREFIX


def jackhmmer_to_fasta(msa):
    """
    Convert the MSA data from jackhmmer to fasta format.

    Args:
        msa (list): Multi-sequence alignment data.

    Returns:
        list: Formatted sequences.
    """
    results_dir = os.path.join("../results",
                               "msas",
                               f"{PREFIX}_hmmer_{len(msa[0])}_seqs.a3m")

    def reformat_sequences(msa):
        formatted_sequences = []
        for idx, sequence in enumerate(msa[0]):
            formatted_sequence = f">sequence_{idx + 1}\n{sequence}\n"
            formatted_sequences.append(formatted_sequence)
        return formatted_sequences

    formatted_sequences = reformat_sequences(msa)

    def save_formatted_sequences_to_file(formatted_sequences):
        with open(results_dir, 'w', encoding="utf-8") as file:
            for sequence in formatted_sequences:
                file.write(sequence)

    save_formatted_sequences_to_file(formatted_sequences)


def run_jackhmmer():
    """
    Execute the jackhmmer command. (Implementation is not provided)
    """
    ### NOT IMPLEMENTED ###
    msa = [['HLDERPLACE', 'PLACEHLDER']]
    jackhmmer_to_fasta(msa)


def build_jackhmmer_msa(sequence):
    """
    Build an MSA for a target sequence using jackhmmer
        (querying UniRef90, small BFD, and mgnify remotely).

    Args:
        sequence (str): Target sequence.
    """
    jackhmmer_installed = os.path.isfile(JACKHMMER_PATH)

    if not jackhmmer_installed:
        print('JACKHMMER not found, not building MSA')
    else:
        run_jackhmmer()
