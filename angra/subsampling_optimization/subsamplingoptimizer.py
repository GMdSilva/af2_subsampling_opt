"""
Defines a class that runs the subsampling parameter optimization pipeline,
    starting at detecting peaks (modes), ranking them based on various criteria,
    measuring the performance of each parameter set, and choosing the top ranked
"""
import os
from typing import Any, Dict, List, Tuple
import json
import numpy as np

from user_settings.config import PREDICTION_ROOT, SYSTEM_NAME
from angra.utilities.plotter import Plotter
from angra.ensemble_analysis.mdanalysisrunner import MDAnalysisRunner
from angra.subsampling_optimization.peakdetector import PeakDetector
from angra.subsampling_optimization.predictionevaluator import PredictionEvaluator
from angra.utilities.utilities import find_common_ranges_with_wiggle,\
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
        self.evaluator = None
        self.final_ranges = None
        self.subsampling_results = None

    def get_final_ranges(self,
                         window_size: int = 10,
                         peak_thresh: int = 2,
                         peak_n: int = 4) -> List[Tuple[int, int]]:
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

        for result in self.subsampling_results:
            res_len = result['residues'].max()
            window_size = int((res_len*0.05))
            detector = PeakDetector(result)
            peak_range = detector.detect_peaks(window_size, peak_thresh, peak_n)
            all_ranges.append(peak_range)
        self.final_ranges = find_common_ranges_with_wiggle(all_ranges)

    def analyze_predictions(self,
                            method: str,
                            trial,
                            selection: str = 'protein and name CA',
                            bulk: bool = True) -> list[dict]:
        """
        Analyze a series of sub-sampled AF2 af2_predictions
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
                                'af2_predictions',
                                f"{self.prefix}"
                                f"_{trial.split(':')[0]}"
                                f"_{trial.split(':')[1]}")
            mdanalyzer = MDAnalysisRunner(path, self.prefix, selection)
            self.subsampling_results = mdanalyzer.process_results(bulk=bulk,
                                        method=method, trial=trial)
        else:
            path = os.path.join(PREDICTION_ROOT,
                                'results',
                                'af2_predictions',
                                '')
            mdanalyzer = MDAnalysisRunner(path, self.prefix, selection)
            self.subsampling_results = mdanalyzer.process_results(bulk=bulk,
                                                             method=method,
                                                                  trial=trial)

    def analyze_parameter_set(self) -> None:
        """
        Compute the variance of differences between adjacent in-between values
        for a bimodal distribution.

        Returns:
            final_ranges: RMSD ranges used for the analysis

        """
        if self.custom_range:
            self.final_ranges = self.selection
        else:
            self.get_final_ranges()

    def make_final_decision(self, path):
        self.evaluator = PredictionEvaluator(self.prefix)
        subsampling_parameters = self.evaluator.final_subsampling_decision(self.final_ranges, path)
        final_results_path = os.path.join(PREDICTION_ROOT,
                                          "results",
                                          "optimization_results",
                                          f"{self.prefix}_optimizer_results.pkl")
        save_to_pickle(final_results_path, subsampling_parameters)
        subsampling_parameters = self._reformat_final_decision(subsampling_parameters)
        return subsampling_parameters

    def plot_rmsd_ranges(self):
        for r_range in self.final_ranges:
            path_results = (os.path.join('results',
                            'misc_data',
                            f"{self.prefix}_"
                            f"{r_range[0]}_"
                            f"{r_range[1]}_"
                            f"modes_all_trials.pkl"))
            mode_data = load_from_pickle(path_results)[1]
            plotter = Plotter(self.prefix, mode_data)
            plotter.plot_kde_with_modes(r_range)

    @staticmethod
    def _reformat_final_decision(final_decision):

        def convert_np_arrays(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()  # Convert ndarray to list
            raise TypeError(f"Object of type '{type(obj).__name__}' is not JSON serializable")

        reformatted = {'System Name': SYSTEM_NAME,
                       'max_seq': final_decision['parameters'].split(':')[0],
                       'extra_seq': final_decision['parameters'].split(':')[1],
                       'Highest Variation Range': final_decision['ranges'][0],
                       'Score': final_decision['scores'][0]}
        path_reports = os.path.join('results',
                                    'reports',
                                    f"{SYSTEM_NAME}_best_subsampling_parameters.json")

        with open(path_reports, 'w') as file:
            json.dump(reformatted, file, default=convert_np_arrays)

        return reformatted

    def plot_rmsf_results(self):
        plotter = Plotter(self.prefix, self.subsampling_results)
        plotter.plot_rmsf_peaks(self.subsampling_results, self.final_ranges)