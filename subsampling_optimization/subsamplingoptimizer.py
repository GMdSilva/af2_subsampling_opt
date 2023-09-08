"""
Defines a class that runs the subsampling parameter optimization pipeline,
    starting at detecting peaks (modes), ranking them based on various criteria,
    measuring the performance of each parameter set, and choosing the top ranked
"""
import os
import re
from typing import Any, Dict, List, Tuple

import matplotlib.pyplot as plt

from ensemble_analysis.mdanalysisrunner import MDAnalysisRunner
from subsampling_optimization.peakdetector import PeakDetector
from subsampling_optimization.predictionevaluator import PredictionEvaluator
from utilities.utilities import find_common_ranges_with_wiggle, load_from_pickle


class SubsamplingOptimizer:
    """
    Class that runs the subsampling parameter optimization pipeline,
        starting at detecting peaks (modes), ranking them based on various criteria,
        measuring the performance of each parameter set, and choosing the top ranked
    """

    def __init__(self, prefix: str,
                 use_custom_range: bool = False,
                 selection='protein and name CA'):
        """
        Args:
        - prefix (str): Usually the name of the target protein.
        - use_custom_range (str): Wether to use custom residue range or not
        - selection (List[Tuple[int, int]]): Custom residue range for analyzes
        """
        self.prefix = prefix
        self.custom_range = use_custom_range
        self.selection = selection

    def plot_peaks(self, result: Dict[str, Any], peak_range: List[Tuple[int, int]]) -> None:
        """
        Plots peaks for the given results.

        Parameters:
        - result: A dictionary containing results and residues.
        - peak_range: A list of tuple ranges indicating peaks.
        """
        colors = ['red', 'blue', 'green', 'yellow', 'purple']
        plt.close()
        plt.plot(result['residues'], result['results'],
                 marker="o",
                 label=result['trial'])
        plt.xlabel('Value (A)', fontsize=18)
        plt.ylabel('Density', fontsize=18)
        plt.title(self.prefix, fontsize=20)
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
                                 f"{self.prefix}_ranges_{result['trial'].split(':')[0]}.png")
        plt.savefig(save_path)

    def get_final_ranges(self, subsampling_results: List[Dict[str, Any]],
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
            detector = PeakDetector(result)
            peak_range = detector.detect_peaks(window_size, peak_thresh, peak_n)
            all_ranges.append(peak_range)
            self.plot_peaks(result, peak_range)

        all_ranges_expanded = find_common_ranges_with_wiggle(all_ranges)

        return all_ranges_expanded

    def analyze_predictions(self,
                            method: str,
                            selection: str = 'protein and name CA',
                            trial: str = "256:512",
                            bulk: bool = True) -> list[dict]:
        """
        Analyze a series of subsampled AF2 predictions
        by calculating MDanalysis observables

        Args:
            method (str): Analysis method (e.g. RMSF, RMSD, etc.).
            selection (str): Residue range to run analysis on.
            trial (str): Subsampling parameter being tested.
            bulk (str): Run analysis per-folder or in bulk.

        Returns:
            subsampling_results (list[dict[str, str]]): list of dicts containing metadata
                about each prediction and the calculated observables
        """
        results_path = os.path.join("results",
                                    "mdanalysis_results",
                                    f"{self.prefix}_{method}_results.pkl")

        file_exists = os.path.isfile(results_path)

        if not bulk:
            path = os.path.join('results',
                                'predictions',
                                f"{self.prefix}"
                                f"_{trial.split(':')[0]}"
                                f"_{trial.split(':')[1]}")
        else:
            path = os.path.join('results',
                                'predictions',
                                '')

        if not file_exists:
            mdanalyzer = MDAnalysisRunner(path, self.prefix, selection)
            subsampling_results = mdanalyzer.process_results(bulk=bulk,
                                                             method=method)
        else:
            subsampling_results = load_from_pickle(results_path)

        return subsampling_results

    def get_optimized_parameters(self, path: str,
                                 subsampling_results: list[dict]) -> dict:
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
                                    f"{self.prefix}_rmsd_results.pkl")
        file_exists = os.path.isfile(results_path)
        if self.custom_range:
            final_ranges = self.selection
        elif file_exists:
            final_ranges = load_from_pickle(results_path)
        else:
            final_ranges = self.get_final_ranges(subsampling_results)
        evaluator = PredictionEvaluator(self.prefix)
        subsampling_parameters = evaluator.final_subsampling_decision(final_ranges, path)

        return subsampling_parameters

