"""
Defines a class that scores each subsampling parameter set based on peak analysis.
"""

from typing import Union, List
import os
import numpy as np

from ensemble_analysis.utils.bulk_ensemble_analysis import bulk_analysis
from subsampling_optimization.peakanalyzer import PeakAnalyzer
from utilities.utilities import tuple_range_to_string, load_from_pickle, save_to_pickle


class PredictionEvaluator:
    """
    Class that scores each subsampling parameter set based on peak analysis
    """
    def __init__(self, prefix: str):
        """
        Args:
        - prefix: Usually, the name of the target protein.
        """
        self.prefix = prefix

    @staticmethod
    def in_between_spacing_metric(data: Union[List[float], np.ndarray],
                                  mode1: float,
                                  mode2: float) -> float:
        """
        Compute the variance of differences between adjacent in-between values
        for a bimodal distribution.

        Args:
            data (Union[List[float], np.ndarray]): The dataset.
            mode1 (float): First mode of the bimodal distribution.
            mode2 (float): Second mode of the bimodal distribution.

        Returns:
            float: Variance of differences between adjacent sorted in-between values.
                   Returns 0.0 if no in-between values are found.

        Notes:
            This function assumes that mode1 < mode2.
        """
        data = np.asarray(data)
        in_between_values = data[(data > mode1) & (data < mode2)]
        if len(in_between_values) == 0:
            return 0.0
        diffs = np.diff(np.sort(in_between_values))
        return np.var(diffs)

    def analyze_distributions(self, distributions: dict, rmsd_range: str, trials: list) -> dict:
        """
        Analyze a set of distributions for modality.

        Parameters:
        - distributions (dict): Dictionary containing distribution results.
        - rmsd_range (str): RMSD range for the analysis.
        - trials (list): List of trial identifiers.

        Returns:
        - dict: Analyzed modalities indexed by trial.
        """
        analyzed_modes = {}
        for i, distribution in enumerate(distributions['results']):
            analyzer = PeakAnalyzer(self.prefix, trials[i])
            bandwidth = analyzer.get_silverman_bandwidth(distribution)
            mode_data = analyzer.analyze_modes(distribution, rmsd_range, bandwidth)
            if mode_data['num_modes'] > 1:
                analyzed_mode_data = analyzer.two_state_analysis(mode_data)
                analyzed_modes[trials[i]] = analyzed_mode_data
        return analyzed_modes

    @staticmethod
    def calculate_metrics(mode_data: dict) -> float:
        """
        Calculate a metric score for given modes.

        Parameters:
        - mode_data (dict): Dictionary containing modality data.

        Returns:
        - float: Calculated metric score.
        """
        intermediate_coverage = PredictionEvaluator.in_between_spacing_metric(
            mode_data['distribution'],
            mode_data['modes'][1],
            mode_data['modes'][0]
        )
        combined_population = mode_data['ground_pop'] + mode_data['alt1_pop']
        distance_between_states = mode_data['dist_ground_alt1']
        score = (intermediate_coverage + distance_between_states) / combined_population
        return score

    def find_largest_score(self,
                           distributions: dict,
                           rmsd_range: str,
                           trials: list) -> tuple:
        """
        Find the trial with the largest score
            after analyzing distributions and performing checks.

        Parameters:
        - distributions (dict): Dictionary containing distribution results.
        - rmsd_range: RMSD range for the analysis.
        - trials (list): List of trial identifiers.

        Returns:
        - tuple: (trial with maximum score, maximum score)
        """
        analyzed_modes = self.analyze_distributions(distributions, rmsd_range, trials)

        scores = {}

        for trial, mode_data in analyzed_modes.items():
            analyzer = PeakAnalyzer(self.prefix, trial)
            if analyzer.mode_sanity_checks(mode_data, rmsd_range):
                score = self.calculate_metrics(mode_data)
                scores[trial] = score

        if scores:
            max_score_trial = max(scores, key=scores.get)
            max_score = scores[max_score_trial]
        else:
            max_score_trial = "failed"
            max_score = 0

        return max_score_trial, max_score

    def choose_subsampling_parameters(self,
                                      loaded_data: dict,
                                      rmsd_range: str,
                                      load_from_disk: bool = False) -> tuple:
        """
        Choose the best subsampling parameters based on the scores.

        Parameters:
        - loaded_data (dict): The data to analyze.
        - rmsd_range (str): The RMSD range to consider.
        - load_from_disk (bool, optional): If data should be loaded from disk. Defaults to False.

        Returns:
        - tuple: Best scoring test and associated scores.
        """
        if load_from_disk:
            loaded_data = load_from_pickle(loaded_data)

        trials = loaded_data['trial']

        best_scoring_test, scores = self.find_largest_score(loaded_data,
                                                            rmsd_range,
                                                            trials)

        return best_scoring_test, scores

    def test_different_peaks(self,
                             final_ranges: list,
                             path: str,
                             save_to_disk: bool = True) -> list:
        """
        Test various peaks and analyze the results.

        Parameters:
        - final_ranges (list): List of peak ranges to test.
        - path (str): Path to the data directory.

        Returns:
        - list: A list of dictionaries containing results for each peak range.
        """
        results = []

        for peak_range in final_ranges:
            selection = f'resid {tuple_range_to_string(peak_range)} and backbone'
            subsampling_results = bulk_analysis(self.prefix, path, selection=selection)

            trials = [d['trial'] for d in subsampling_results if 'trial' in d]
            values = [d['results'] for d in subsampling_results if 'results' in d]

            result_dict_partial = {
                'results': values,
                'trial': trials,
            }

            best_scoring, score = self.choose_subsampling_parameters(result_dict_partial,
                                                                     peak_range)

            result_dict = {
                'range': peak_range,
                'results': result_dict_partial,
                'chosen_parameters': best_scoring,
                'score': score,
            }
            results.append(result_dict)

        if save_to_disk:
            results_path = os.path.join('results',
                                        'mdanalysis_results',
                                        f"{self.prefix}_peak_selection_results.pkl")
            save_to_pickle(results_path, results)

        return results

    @staticmethod
    def create_parameters_dict(parameters: str, ranges: list, scores: list) -> dict:
        """Create a dictionary containing subsampling parameters."""
        subsampling_parameters = parameters.split(':')
        sorted_indexes = np.argsort(scores)[::-1]

        return {
            'parameters': subsampling_parameters[1],
            'ranges': np.array(ranges)[sorted_indexes],
            'scores': np.array(scores)[sorted_indexes]
        }

    def final_subsampling_decision(self, final_ranges: list, path: str) -> dict:
        """
        Make a decision on the best subsampling parameters based on results from various peaks.

        Parameters:
        - final_ranges (list): List of peak ranges.
        - path (str): Path to the data directory.

        Returns:
        - dict: Dictionary containing the chosen parameters, their ranges, and scores.
        """
        results_path = os.path.join("results",
                                    "mdanalysis_results",
                                    f"{self.prefix}_peak_selection_results.pkl")

        file_exists = os.path.isfile(results_path)
        if not file_exists:
            results = self.test_different_peaks(final_ranges, path)
        else:
            results = load_from_pickle(results_path)

        all_same_results = all(x['chosen_parameters']
                               == results[0]['chosen_parameters']
                               for x in results)
        scores = [d['score'] for d in results if 'score' in d]

        if all_same_results and results[0]['chosen_parameters'] == 'failed':
            print('Failed :(')
            return {}

        if all_same_results:
            print(f'Unanimous Success! '
                  f'Choosing {results[0]["chosen_parameters"].split(":")[0]}'
                  f' as subsampling parameters')
            return self.create_parameters_dict(results[0]['chosen_parameters'],
                                               final_ranges,
                                               scores)

        if not all_same_results:
            index = scores.index(max(scores))
            print(f'Will use {results[index]["chosen_parameters"].split(":")[0]}'
                  f' for analyzes')
            return self.create_parameters_dict(results[index]['chosen_parameters'],
                                               final_ranges,
                                               scores)

        print('Failed :(')
        return {}
