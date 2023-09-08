from typing import Dict
import re
import os

from subsampling_optimization.subsamplingoptimizer import SubsamplingOptimizer
from subsampling_optimization.predictionevaluator import PredictionEvaluator
from utilities.utilities import load_from_pickle, save_to_pickle

class MutationAnalyzer(SubsamplingOptimizer):
    def __init__(self, prefix, selection):
        self.prefix = prefix
        self.selection = selection

    def compare_conditions(self, trials, old_prefix: str = 'abl_wt') -> Dict[str, float]:
        """
        Compares conditions based on a given path and trials.

        Parameters:
        - trials: Another parameter (provide a meaningful description)
        - old_prefix (str, optional): Prefix for older data. Defaults to 'abl_wt'.

        Returns:
        - differences (Dict[str:float]): Changes to relative state populations per mutation.
        """

        residue_range = self._extract_residue_range()
        subsampling_results = self.analyze_predictions('rmsd', self.selection, bulk=False)
        differences = self.contrast_differences(subsampling_results[0],
                                                     residue_range,
                                                     trials,
                                                     old_prefix)
        return differences

    def contrast_differences(self,
                             distribution: dict,
                             rmsd_range: str,
                             trials: str,
                             old_prefix: str = 'abl_wt') -> Dict[str, float]:
        """
        Analyze and print the contrast differences between a given distribution and a control.

        Parameters:
        - distribution (dict): The distribution to be analyzed.
        - rmsd_range (str): Residue range used to calculate the RMSD.
        - trials (str): Trials to be processed.
        - old_prefix (str, optional): Prefix for the control naming. Defaults to 'abl_wt'.

        Returns:
        - diff_dict Dict[str: float]: Dict containing the differences
            between analyzed states from test and control predictions.
        """

        # Analyze the distribution
        evaluator = PredictionEvaluator(self.prefix)
        analyzed_modes = evaluator.analyze_distribution(distribution, rmsd_range[0], trials)

        # Extract necessary components from the first trial
        trial_part1, trial_part2 = trials[0][0].split(':')

        # Construct control name
        control_name = (f"accepted_{trial_part1}_{trial_part2}_"
                        f"{old_prefix}_"
                        f"range_{rmsd_range[0][0]}_{rmsd_range[0][1]}")

        # Load control data
        control_path = os.path.join("results", "peaks", control_name)
        control = load_from_pickle(control_path)

        # Extract values from the test and control datasets
        ground_pop_test = analyzed_modes[trials[0][0]]['ground_pop']
        alt1_pop_test = analyzed_modes[trials[0][0]]['alt1_pop']
        in_between_pop_test = 1 - ground_pop_test - alt1_pop_test

        ground_pop_control = control['ground_pop']
        alt1_pop_control = control['alt1_pop']
        in_between_pop_control = 1 - ground_pop_control - alt1_pop_control

        # Calculate the differences
        ground_pop_diff = ground_pop_test - ground_pop_control
        alt1_pop_diff = alt1_pop_test - alt1_pop_control
        in_between_diff = in_between_pop_test - in_between_pop_control

        # Construct the result dictionary
        diff_dict = {
            self.prefix: {
                'ground_pop_diff': ground_pop_diff*100,
                'alt1_pop_diff': alt1_pop_diff*100,
                'in_between_diff': in_between_diff*100,
                'ground_pop_test': ground_pop_test*100,
                'alt1_pop_test': alt1_pop_test*100,
                'in_between_pop_test': in_between_pop_test*100,
                'ground_pop_control': ground_pop_control*100,
                'alt1_pop_control': alt1_pop_control*100,
                'in_between_pop_control': in_between_pop_control*100
            }
        }
        return diff_dict

    def _extract_residue_range(self):
        """
        Extracts the residue range from the selection attribute.

        Returns:
        list of tuple: A list containing a single tuple
            that denotes the start and end of the residue range.
        """

        match = re.search(r'(\d+):(\d+)', self.selection)
        if match:
            start, end = map(int, match.groups())
            return [(start, end)]
        return []

    def plot_results(self, list_of_dicts, x_values):
        import matplotlib.pyplot as plt

        # Extract prefixes and ground_pop_diff values
        prefixes = [list(d.keys())[0] for d in list_of_dicts]
        ground_pop_diffs = [d[prefix][x_values] for d in list_of_dicts for prefix in d]

        # Plotting
        fig, ax = plt.subplots()
        ax.barh(prefixes, ground_pop_diffs)
        ax.set_xlabel(x_values)
        ax.set_ylabel('Prefix')
        ax.set_title('')

        plt.tight_layout()
        save_path = os.path.join('results',
                                 'plots',
                                 f"{self.prefix}_{x_values}.png")
        plt.savefig(save_path)
