"""
Functions for building MSA for subsampled AF2 predictions
"""

import os
from user_settings.config import COLABFOLDBATCH_PATH, PREFIX


def estimate_starting_parameters(msa: str) -> int:
    """
    Estimate starting parameters based on MSA depth.

    Args:
    - msa (str): Input multi-sequence alignment.

    Returns:
    - int: Starting exponent for sequence count.
    """
    msa_depth = int(msa.split("_")[-2])
    if msa_depth <= 2000:
        return 2
    if 2000 < msa_depth <= 20000:
        return 3
    return 4


def build_parameter_range(trials: int, msa: str) -> list:
    """
    Construct a range of parameters for AF2 based on trials and MSA.

    Args:
    - trials (int): Number of trials to perform.
    - msa (str): Input multi-sequence alignment.

    Returns:
    - list: List of parameter sets.
    """
    starting_max_seq_exp = estimate_starting_parameters(msa)
    parameters = [
        f"{2 ** exp}:{2 ** (exp + 1)}"
        for exp in range(starting_max_seq_exp, trials + starting_max_seq_exp)
    ]
    return parameters


def run_af2(parameters: list, msa: str, seeds: int, kind: str):
    """
    Execute the AF2 process with the given parameters.

    Args:
    - parameters (list): List of parameter sets.
    - msa (str): Input multi-sequence alignment.
    - seeds (int): Number of seeds.
    - kind (str): Type for prediction.
    """
    result_path = os.path.join('../results', 'predictions')
    for p_set in parameters:
        max_seq, extra_seq = p_set.split(':')
        command = (
            f"colabfold_batch --num-seeds {seeds} --max-msa {p_set} --use-dropout "
            f"{msa} {os.path.join(result_path, f'{PREFIX}_{kind}_{max_seq}_{extra_seq}')}"
        )
        os.system(command)


def make_predictions(msa: str, trials: int, seeds: int, kind: str):
    """
    Generate predictions using the given MSA, trials, seeds, and type.

    Args:
    - msa (str): Input multi-sequence alignment.
    - trials (int): Number of trials to perform.
    - seeds (int): Number of seeds.
    - kind (str): Type for prediction.
    """
    if os.path.isfile(COLABFOLDBATCH_PATH):
        print('colabfoldbatch not found, not making predictions')
        return

    parameters = build_parameter_range(trials, msa)
    run_af2(parameters, msa, seeds, kind)


def run_subsampled_af2(msa: str, trials: int = 5, seeds: int = 32, kind: str = 'wt'):
    """
    Entry point function to generate subsampled AF2 predictions.

    Args:
    - msa (str): Input multi-sequence alignment.
    - trials (int, optional): Number of trials to perform. Default is 5.
    - seeds (int, optional): Number of seeds. Default is 32.
    - kind (str, optional): Type for prediction. Default is 'wt'.
    """
    make_predictions(msa, trials, seeds, kind)
