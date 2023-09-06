"""Finds the residue ranges for calculating distributions for scoring subsampling parameters"""

import os.path
from typing import List, Dict, Tuple, Any

import matplotlib.pyplot as plt

from subsampling_opt.peak_calling import detect_peaks
from utilities.utilities import find_common_ranges_with_wiggle
from user_settings.config import PREFIX


def plot_peaks(result: Dict[str, Any], peak_range: List[Tuple[int, int]]) -> None:
    """
    Plots peaks for the given results.

    Parameters:
    - result: A dictionary containing results and residues.
    - peak_range: A list of tuple ranges indicating peaks.

    Returns:
    - None
    """
    colors = ['red', 'blue', 'green', 'yellow', 'purple']
    plt.close()
    plt.plot(result['residues'], result['results'],
             marker="o",
             label=result['trial'])
    plt.xlabel('Value (A)', fontsize=18)
    plt.ylabel('Density', fontsize=18)
    plt.title(f'{str(PREFIX)}', fontsize=20)
    plt.legend()

    for i, (start, end) in enumerate(peak_range):
        color = colors[i % len(colors)]
        plt.plot(result['residues'][start:end],
                 result['results'][start:end],
                 marker='o',
                 color=color, linestyle='-')
    plt.tight_layout()
    save_path = os.path.join('results',
                             'plots',
                             f"{PREFIX}_ranges_{result['trial'].split(':')[0]}.png")
    plt.savefig(save_path)


def get_final_ranges(subsampling_results: List[Dict[str, Any]],
                     window_size: int = 5,
                     peak_thresh: int = 5,
                     peak_n: int = 3) -> List[Tuple[int, int]]:
    """
    Gets final ranges from subsampling results.

    Parameters:
    - subsampling_results: List of dictionaries containing results, residues, and trials.
    - window_size: Window size for peak detection.
    - peak_thresh: Threshold for peak detection.
    - peak_n: Number of peaks for detection.

    Returns:
    - List of tuple ranges indicating the final expanded ranges.
    """
    all_ranges = []

    for result in subsampling_results:
        peak_range = detect_peaks(result, window_size, peak_thresh, peak_n)
        all_ranges.append(peak_range)
        plot_peaks(result, peak_range)

    all_ranges_expanded = find_common_ranges_with_wiggle(all_ranges)

    return all_ranges_expanded
