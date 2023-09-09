"""
Defines a class that runs the subsampling parameter optimization pipeline,
    starting at detecting peaks (modes), ranking them based on various criteria,
    measuring the performance of each parameter set, and choosing the top ranked
"""
import os
from typing import Any, Dict, List, Tuple

from user_settings.config import PREDICTION_ROOT
from src.utilities.plotter import Plotter
from src.ensemble_analysis.mdanalysisrunner import MDAnalysisRunner
from src.subsampling_optimization.peakdetector import PeakDetector
from src.subsampling_optimization.predictionevaluator import PredictionEvaluator
from src.utilities.utilities import find_common_ranges_with_wiggle,\
    load_from_pickle,\
    save_to_pickle


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
        - use_custom_range (str): To use custom residue range or not
        - selection (List[Tuple[int, int]]): Custom residue range for analyzes
        """
        self.prefix = prefix
        self.custom_range = use_custom_range
        self.selection = selection

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
            plotter = Plotter(self.prefix, result)
            plotter.plot_rmsf_peaks(peak_range)

        all_ranges_expanded = find_common_ranges_with_wiggle(all_ranges)

        return all_ranges_expanded

    def analyze_predictions(self,
                            method: str,
                            selection: str = 'protein and name CA',
                            trial: str = "256:512",
                            bulk: bool = True) -> list[dict]:
        """
        Analyze a series of sub-sampled AF2 predictions
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
        if not bulk:
            path = os.path.join(PREDICTION_ROOT,
                                'results',
                                'optimizer_results',
                                f"{self.prefix}"
                                f"_{trial.split(':')[0]}_"
                                f"_{trial.split(':')[1]}.pkl")
            file_exists = os.path.isfile(path)

            if file_exists:
                subsampling_results = load_from_pickle(path)
                return subsampling_results
            path = os.path.join(PREDICTION_ROOT,
                                'results',
                                'predictions',
                                f"{self.prefix}"
                                f"_{trial.split(':')[0]}_"
                                f"_{trial.split(':')[1]}")
            mdanalyzer = MDAnalysisRunner(path, self.prefix, selection)
            subsampling_results = mdanalyzer.process_results(bulk=bulk,
                                                             method=method)
        else:
            path = os.path.join(PREDICTION_ROOT,
                                "results",
                                "optimizer_results",
                                f"{self.prefix}_{method}_all_results.pkl")

            file_exists = os.path.isfile(path)
            if file_exists:
                subsampling_results = load_from_pickle(path)
                return subsampling_results
            path = os.path.join(PREDICTION_ROOT,
                                'results',
                                'predictions',
                                '')

            mdanalyzer = MDAnalysisRunner(path, self.prefix, selection)
            subsampling_results = mdanalyzer.process_results(bulk=bulk,
                                                             method=method)

        return subsampling_results

    def get_optimized_parameters(self, path: str,
                                 subsampling_results: list[dict]) -> None:
        """
        Compute the variance of differences between adjacent in-between values
        for a bimodal distribution.

        Args:
            path (str): Path to directory containing AF2 subsampling directories
            subsampling_results (str): metadata and prediction observables
                calculated with MDanalysis

        Returns:
            None: saves results to disk.

        """
        results_path = os.path.join(PREDICTION_ROOT,
                                    "results",
                                    "optimizer_results",
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
        final_results_path = os.path.join(PREDICTION_ROOT,
                                          "results",
                                          "optimizer_results",
                                          f"{self.prefix}_optimizer_results.pkl")
        save_to_pickle(final_results_path, subsampling_parameters)
