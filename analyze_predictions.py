"""Main functions for analyzing the output of subsampled AF2"""

import os

from ensemble_analysis.bulk_ensemble_analysis import bulk_analysis
from subsampling_opt.parameter_scoring import final_subsampling_decision
from subsampling_opt.find_ranges import get_final_ranges
from utilities.utilities import load_from_pickle
from user_settings.config import PREFIX


def analyze_predictions(path: str,
                        method: str) -> list[dict]:
    """
    Analyze a series of subsampled AF2 predictions
    by calculating MDanalysis observables

    Args:
        path (str): Path to directory containing AF2 subsampling directories
        method (str): Analysis method (e.g. RMSF, RMSD, etc.)

    Returns:
        subsampling_results (list[dict[str, str]]): list of dicts containing metadata
            about each prediction and the calculated observables
    """

    results_path = os.path.join("results",
                                "mdanalysis_results",
                                f"{PREFIX}_{method}_results.pkl")

    file_exists = os.path.isfile(results_path)

    if not file_exists:
        subsampling_results = bulk_analysis(path, method)
    else:
        subsampling_results = load_from_pickle(results_path)

    return subsampling_results


def get_optimized_parameters(path: str, subsampling_results: list[dict]) -> dict:
    """
    Compute the variance of differences between adjacent in-between values
    for a bimodal distribution.

    Args:
        path (str): Path to directory containing AF2 subsampling directories
        subsampling_results (str): metadata and prediction observables
            calculated with MDanalysis

    Returns:
        subsampling_parameters (dict): dict containing the optimized
            subsampling parameters and metadata about the prediction

    """
    results_path = os.path.join("results",
                                "mdanalysis_results",
                                f"{PREFIX}_rmsd_results.pkl")
    file_exists = os.path.isfile(results_path)

    if not file_exists:
        final_ranges = get_final_ranges(subsampling_results)
    else:
        final_ranges = load_from_pickle(results_path)
    subsampling_parameters = final_subsampling_decision(final_ranges, path)

    return subsampling_parameters
