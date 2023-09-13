"""
Defines class for analyzing mutations and measuring accuracy
    at predicting the effects of mutations in the ground or alternative states
"""

from typing import Any, Dict, List
import os

from src.subsampling_optimization.subsamplingoptimizer import SubsamplingOptimizer
from src.subsampling_optimization.predictionevaluator import PredictionEvaluator
from src.utilities.utilities import load_from_pickle, save_to_pickle
from src.utilities.plotter import Plotter
from src.utilities.utilities import range_to_string
from user_settings.new_config import load_config
from user_settings import config


class MutationAnalyzer(SubsamplingOptimizer):
    """
    Class for analyzing mutations and measuring accuracy
        at predicting the effects of mutations in the ground or alternative states
    """
    def __init__(self, prefix: str):
        """
        Compares conditions based on a given path and trials.

        Args:
        - prefix (str): Usually the name of the protein being studied.
        """
        self.prefix = prefix
        self.selection: str = 'protein and name CA'
        self.wildtype_results: dict = None

    def compare_conditions(self, trial: str) -> Dict[str, float]:
        """
        Compares results of different prediction conditions
            (usually wild-type vs. mutants)

        Parameters:
        - trial: Subsampling parameters used for predicting both conditions.

        Returns:
        - differences (Dict[str:float]): Changes to relative state populations
            per mutation.
        """
        residue_range = self.selection
        selection = f'resid {range_to_string(self.selection)} and name CA'
        subsampling_results = self.analyze_predictions('rmsd', selection, bulk=False)
        differences = self.contrast_differences(subsampling_results,
                                                residue_range,
                                                trial)
        return differences

    def get_control_values(self, trial: str,
                           residue_range: str) -> dict:
        """
        Retrieve wild-type values based on the trials and rmsd_range.

        Parameters:
        - trial (str): Subsampling parameters used for predicting both conditions.
        - residue_range (str): Residue range used to calculate the RMSD.
            Defaults to 'abl_wt'.
        """

        # Extract necessary components from the first trial
        trial_part1, trial_part2 = trial[0][0].split(':')

        # Construct control name
        wildtype_filename = (f"{config.SYSTEM_NAME}_accepted_{trial_part1}_{trial_part2}_"
                             f"range_{residue_range[0]}_{residue_range[1]}.pkl")

        # Load control data
        wildtype_path = os.path.join(config.PREDICTION_ROOT,
                                     "results",
                                     "optimization_results",
                                     "peaks",
                                     wildtype_filename)
        self.wildtype_results = load_from_pickle(wildtype_path)

    def contrast_differences(self,
                             distribution: dict,
                             residue_range: str,
                             trial: str) -> dict:
        """
        Compares results of different prediction conditions
            (usually wild-type vs. mutants)

        Parameters:
        - distribution (dict): Subsampling parameters used for predicting both conditions.
        - residue_range (str): Prefix for older data. Defaults to 'abl_wt'.
        - trial (str): Subsampling parameters used for predicting both conditions.

        Returns:
        - differences (Dict[str:float]): Changes to relative state populations
            per mutation.
        """
        # Analyze the distribution
        evaluator = PredictionEvaluator(self.prefix)
        analyzed_modes = evaluator.analyze_distributions(distribution,
                                                         residue_range,
                                                         trial)

        self.get_control_values(trial, residue_range)

        # Extract values from the test and control datasets
        ground_pop_test = analyzed_modes[trial[0][0]]['ground_pop']
        alt1_pop_test = analyzed_modes[trial[0][0]]['alt1_pop']
        in_between_pop_test = 1 - ground_pop_test - alt1_pop_test

        ground_pop_control = self.wildtype_results['ground_pop']
        alt1_pop_control = self.wildtype_results['alt1_pop']
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

    def load_or_generate_mut_results(self,
                                     results_path: str,
                                     optimization_results: dict)\
            -> List[Dict[str, Any]]:
        """
        Loads or generates mutation vs. wildtype comparison
            for a list of mutants.

        Args:
            results_path (str): Path where to save/read results
            optimization_results (Dict[str, Any]): The results from
                the subsampling optimization.

        Returns:
            all_mut_results: List[Dict[str, Any]]:
                Results from the mutation vs. wt comparison.
            mut_data: Dict[str, Any]]: sorted mutation metadata
        """
        all_mut_results = []
        file_exists = os.path.isfile(results_path)

        mut_data = config.MUTANT_DATA
        mut_data = dict(sorted(mut_data.items(), key=lambda item: item[1]["rank"]))

        if file_exists:
            all_mut_results = load_from_pickle(results_path)
            self.plot_mut_results(all_mut_results, mut_data)
            return all_mut_results, mut_data

        self.selection = optimization_results['ranges'][0]

        for mut in mut_data:
            self.prefix = mut
            trial = [[(optimization_results['parameters'])]]
            mut_results = self.compare_conditions(trial=trial)
            all_mut_results.append(mut_results)

        self.plot_mut_results(all_mut_results, mut_data)
        save_to_pickle(results_path, all_mut_results)
        return all_mut_results, mut_data

    @staticmethod
    def measure_accuracy(mut_results: dict, mut_data: dict)\
            -> List:
        """
        Evaluates prediction results and measures accuracy
            for predicting the effects of mutations
            on relative state populations.

        Args:
            mut_results (dict): Results from the mutant/wt comparisons
            mut_data (dict): Mutant metadata

        Returns:
            List[Dict, int, int]: Results for each prediction,
                plus % accuracy for each state.
        """
        def calculate_accuracy(diff_key: str,
                               effect_key: str,
                               data: dict,
                               predictions_dict: dict) -> float:
            """
            Calculates accuracy for total set of af2_predictions

            Args:
                diff_key (str): Key to calculate accuracy on.
                effect_key (str): Key that defines the effect of the mutation.
                data (dict): Dict containing mutation analysis results.
                predictions_dict (dict): Dict containing mutation metadata.

            Returns:
                accuracy (float): Accuracy % for given prediction set.
            """
            accuracy = 0
            right_predictions, wrong_predictions = {}, {}
            for entry in data:
                for key, _ in entry.items():
                    if mut_data[key]['effect'][effect_key] == 'ref':
                        continue

                    diff = entry[key][diff_key]
                    effect = mut_data[key]['effect'][effect_key]

                    if (diff > 0 and effect == '+') or (diff < 0 and effect == '-'):
                        accuracy += 1
                        right_predictions[key] = f"Right by {diff}%"
                    else:
                        token = '+' if diff > 0 else '-'
                        wrong_predictions[key] = f"Should be {effect}, is {token}, wrong by {diff}%"

            predictions_dict[effect_key + ' right'] = right_predictions
            predictions_dict[effect_key + ' wrong'] = wrong_predictions
            return accuracy

        total_measurements = len(mut_results) - 1
        predictions = {}

        ground_pop_accuracy = calculate_accuracy('ground_pop_diff',
                                                 'ground_pop',
                                                 mut_results,
                                                 predictions)
        alt1_pop_accuracy = calculate_accuracy('alt1_pop_diff',
                                               'alt1_pop',
                                               mut_results,
                                               predictions)

        return(predictions,
               ground_pop_accuracy/total_measurements*100,
               alt1_pop_accuracy/total_measurements*100)

    def plot_mut_results(self,
                         all_mut_results: List[Dict[str, Any]],
                         mut_data: dict) -> None:
        """
        Plots mutation results for a given prefix and range.

        Args:
            all_mut_results (List[Dict[str, Any]]): A list of mutation results.
            mut_data (dict): Metadata for mutations.

        Returns:
            None
        """
        labels = {'ground_pop_diff': "Ground State Pop. Δ (%)",
                  'alt1_pop_diff': "Alt1 State Pop. Δ (%)",
                  'ground_pop_test': "Ground State Pop. (%)",
                  'alt1_pop_test': "Alt1 State Pop. (%)"}
        plotter = Plotter(self.prefix, all_mut_results)
        for key, value in labels.items():
            pair = {key: value}
            plotter.plot_mut_analysis(pair, mut_data)
