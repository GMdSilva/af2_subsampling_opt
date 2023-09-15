"""
Defines a class that finds and saves representative structures for each selected peak.
"""

from typing import Dict, Any, Tuple, List, Union
import os
import glob
import shutil
import numpy as np

from angra.utilities.utilities import load_from_pickle
from angra.ensemble_analysis.mdanalysisrunner import MDAnalysisRunner
from angra.utilities.plotter import Plotter
from user_settings import config


class StateFinder:
    """
    Class for analyzing mutations and measuring accuracy at predicting the effects
    of mutations in the ground or alternative states.
    """
    def __init__(self, prefix: str):
        """
        Initializes the MutationAnalyzer class.

        Args:
            prefix (str): Name of the protein being studied.
        """
        self.prefix = prefix
        file = load_from_pickle(os.path.join(config.PREDICTION_ROOT,
                                                        'results',
                                                        "optimization_results",
                                                        f"{self.prefix}_accepted_trials.pkl"))
        self.selection: str = 'protein and name CA'
        self.trial: str = None

        d = {}
        for sublist in file:
            key = sublist[0]
            if key not in d or sublist[2] > d[key][2]:
                d[key] = sublist

        result = list(d.values())

        self.all_trials = [sublist[0] for sublist in result]
        all_scores = [sublist[2] for sublist in result]

        self.annot = []
        for score in all_scores:
            self.annot.append(f"Score: {score}")


    def get_variables(self, optimization_results: Dict[str, Any]) -> None:
        """
        Sets trial and range based on given results path and optimization results.

        Args:
            optimization_results (dict): Dictionary containing optimization results.
        """
        self.selection = f"{optimization_results['ranges'][0][0]}_" \
                         f"{optimization_results['ranges'][0][1]}"
        self.trial = optimization_results['parameters'].replace(":", "_")

    def get_reps(self, optimization_results: Dict[str, Any]) -> dict:
        """
        Finds the representative structures for each peak.
        """
        self.get_variables(optimization_results)
        peak_calling_path = os.path.join(
            config.PREDICTION_ROOT,
            'results',
            'optimization_results',
            'peaks',
            f"{self.prefix}_accepted_{self.trial}"
            f"_range_{self.selection}.pkl")
        peak_calling_results = load_from_pickle(peak_calling_path)
        rmsd_path = os.path.join(
            config.PREDICTION_ROOT,
            'results',
            'misc_data',
            f"{self.prefix}_rmsd_{self.selection}_noref_{self.trial}.pkl")
        rmsd_results = load_from_pickle(rmsd_path)

        # Finding index for ground
        index_ground = StateFinder._find_closest_index(
            rmsd_results['results'], peak_calling_results['modes'][-1])
        prediction_ground = rmsd_results['residues'][index_ground]

        # Finding index for alt
        index_alt = StateFinder._find_closest_index(
            rmsd_results['results'], peak_calling_results['modes'][-2])
        prediction_alt1 = rmsd_results['residues'][index_alt]
        return {'ground_ref_index': int(prediction_ground),
                'alt1_ref_index': int(prediction_alt1)}

    @staticmethod
    def _find_closest_index(results: np.ndarray, target_value: float) -> int:
        """
        Finds the index of the value in results that is closest to target_value.

        Args:
            results (np.ndarray): Array of values.
            target_value (float): Target value to find the closest match.

        Returns:
            int: Index of the closest value.
        """
        return np.argmin(np.abs(results - target_value))

    def get_refs_and_compare(self,
                             optimization_results: Dict[str, int]) -> Tuple[List[str], List[str]]:
        """
        Finds PDB files based on optimization results.

        Args:
            optimization_results (Dict[str, int]): A dictionary containing indexes
                                                  for ground and alt1 references.

        Returns:
            Tuple[List[str], List[str]]: A tuple containing two lists.
                                         The first list contains matching paths for ground,
                                         and the second list contains matching paths for alt1.
        """
        indexes = self.get_reps(optimization_results)

        ground_ref_path = self._construct_path(indexes['ground_ref_index'])
        alt1_ref_path = self._construct_path(indexes['alt1_ref_index'])

        ground_files = glob.glob(ground_ref_path)

        alt1_files = glob.glob(alt1_ref_path)

        labeled_files = list(zip(ground_files,
                                 ["ground"] * len(ground_files))) + \
                        list(zip(alt1_files, ["alt1"] * len(alt1_files)))

        rmsds = {}
        for file, label in labeled_files:
            self._save_rep_states(file, label)
            rmsd = self.get_trials_rmsd_vs_refs(file, label)
            rmsds[label] = rmsd

        self._plot_results_trials(rmsds)
        return rmsds

    def _construct_path(self, index: int) -> str:
        """
        Constructs a file path based on a given index.

        Args:
            index (int): The index used in the file name.

        Returns:
            str: The constructed file path.
        """
        index += 1
        formatted_index = str(index).zfill(3)
        return os.path.join('results', 'af2_predictions',
                            f"{self.prefix}_{self.trial}",
                            f"{self.prefix}_unrelaxed_rank_{formatted_index}*.pdb")

    def _save_rep_states(self, file: str, label: str) -> None:
        """
        Save representative states.

        Parameters:
        - file (str): Path to the file.
        - label (str): Label for the state.

        Returns:
        - None
        """
        dest_folder = os.path.join('results', 'references',
                                   f"{self.prefix}_{self.trial}_{self.selection}_{label}_ref.pdb")
        shutil.copy2(file, dest_folder)

    def get_rmsd_vs_refs(self, ref: str,
                         label: str) -> Dict[str, Union[str, List[str]]]:
        """
        Get RMSD against reference.

        Parameters:
        - file (str): Path to the file.
        - label (str): Label for the state.

        Returns:
        - Dict[str, Union[str, List[str]]]: Dictionary of results.
        """
        prefixes = self._get_prefixes()
        path = os.path.join('results', 'af2_predictions')
        results = {}
        trial = self.trial.replace("_", ":")
        selection = f"resid {self.selection.replace('_', ':')} and name CA"

        for prefix in prefixes:
            file_exists = self._check_if_file_exists(prefix, label, self.trial)
            if file_exists:
                res = file_exists
            else:
                runner = MDAnalysisRunner(path,
                                          prefix,
                                          selection,
                                          ref_name=ref)
                res = runner.process_results(bulk=False,
                                             trial=trial,
                                             method='rmsd_ref')
            results[prefix] = res
        return results

    def get_trials_rmsd_vs_refs(self, file: str, label: str,) -> Dict[str, Union[str, List[str]]]:
        """
        Get RMSD against reference.

        Parameters:
        - file (str): Path to the file.
        - label (str): Label for the state.

        Returns:
        - Dict[str, Union[str, List[str]]]: Dictionary of results.
        """
        ref = file
        path = os.path.join('results', 'af2_predictions')
        results = {}
        trials = []

        selection = f"resid {self.selection.replace('_', ':')} and name CA"
        for trial in self.all_trials:
            trial_r = trial.replace("_", ":")
            file_exists = self._check_if_file_exists(self.prefix, label, trial)
            if file_exists:
                res = file_exists
            else:
                runner = MDAnalysisRunner(path,
                                          self.prefix,
                                          selection,
                                          ref_name=ref)
                res = runner.process_results(bulk=False,
                                             trial=trial_r,
                                             method='rmsd_ref')
            results[trial] = res
        return results

    def _build_results_dict_trials(self, rmsd_dict: Dict) -> Dict:
        """
        Build results dictionary from RMSD dictionary.

        Parameters:
        - rmsd_dict (Dict): RMSD dictionary.

        Returns:
        - Dict: Formatted results dictionary.
        """
        results_tuples, ranks = [], []
        labels = StateFinder._get_mutation_names()
        for trial in self.all_trials:
            result_tuple = (rmsd_dict['ground'][trial][0]['results'],
                            rmsd_dict['alt1'][trial][0]['results'])
            rank = [int(value) for value in rmsd_dict['ground'][trial][0]['residues']]
            ranks.append(rank)
            results_tuples.append(result_tuple)

        return {'results_tuples': results_tuples,
                'results_ranks': ranks,
                'results_labels': self.all_trials}

    def _build_results_dict(self, rmsd_dict: Dict) -> Dict:
        """
        Build results dictionary from RMSD dictionary.

        Parameters:
        - rmsd_dict (Dict): RMSD dictionary.

        Returns:
        - Dict: Formatted results dictionary.
        """
        labels = StateFinder._get_mutation_names()
        prefixes = StateFinder._get_prefixes()
        effects = StateFinder._get_effects()

        results_tuples, ranks = [], []

        for prefix in prefixes:
            result_tuple = (rmsd_dict['ground'][prefix][0]['results'],
                            rmsd_dict['alt1'][prefix][0]['results'])
            rank = [int(value) for value in rmsd_dict['ground'][prefix][0]['residues']]
            ranks.append(rank)
            results_tuples.append(result_tuple)

        return {'results_tuples': results_tuples,
                'results_ranks': ranks,
                'results_labels': labels,
                'results_effects': effects}

    def _check_if_file_exists(self, prefix: str, label: str, trial: str) -> Union[bool, Dict]:
        """
        Check if a specific file exists.

        Parameters:
        - prefix (str): Prefix of the file.
        - label (str): Label for the file.

        Returns:
        - Union[bool, Dict]: F if file doesn't exist, else returns the file data.
        """
        file_path = os.path.join(config.PREDICTION_ROOT, 'results', 'misc_data',
                                 f"{prefix}_rmsd_ref_"
                                 f"{trial}_{self.selection}_{label}_"
                                 f"results.pkl")
        if os.path.exists(file_path):
            return load_from_pickle(file_path)
        return False

    @staticmethod
    def _sort_mut_list() -> Dict:
        """
        Sort mutation list from configuration.

        Returns:
        - Dict: Sorted mutation list.
        """
        mut_list = config.MUTANT_DATA
        return dict(sorted(mut_list.items(), key=lambda item: item[1]["rank"]))

    @staticmethod
    def _get_prefixes() -> List[str]:
        """
        Extract prefixes from the sorted mutation list.

        Returns:
        - List[str]: List of prefixes.
        """
        return list(StateFinder._sort_mut_list().keys())

    @staticmethod
    def _get_effects() -> List[str]:
        """
        Extract effects from the sorted mutation list.

        Returns:
        - List[str]: List of effect labels or None if no effect is found.
        """
        muts_list = StateFinder._sort_mut_list()
        prefixes = StateFinder._get_prefixes()

        effects = []
        for prefix in prefixes:
            effect = muts_list[prefix].get('effect')
            if not effect:
                effects.append(f"Ground: {effect['ground_pop']} Alt1: {effect['alt1_pop']}")
            else:
                effects.append('')
        return effects

    @staticmethod
    def _get_mutation_names() -> List[str]:
        """
        Extract mutation names from the sorted mutation list.

        Returns:
        - List[str]: List of mutation names.
        """
        muts_list = StateFinder._sort_mut_list()
        return [muts_list[prefix]['label'] for prefix in StateFinder._get_prefixes()]

    def _plot_results(self, rmsds) -> None:
        """
        Plot results using RMSD data.

        Parameters:
        - rmsds (Dict): RMSD data.

        Returns:
        - None
        """
        results_dict = self._build_results_dict(rmsds)
        plotter = Plotter(self.prefix, results_dict['results_tuples'])
        plotter.plot_multiple_scatter(data_list=results_dict['results_tuples'],
                                      labels=results_dict['results_labels'],
                                      annotations=results_dict['results_effects'],
                                      colors=results_dict['results_ranks'],
                                      filename='ground_vs_alt1_plddt_ranked_muts',
                                      colorbar_label='Ranking')

    def _plot_results_trials(self, rmsds: Dict) -> None:
        """
        Plot results using RMSD data.

        Parameters:
        - rmsds (Dict): RMSD data.

        Returns:
        - None
        """
        results_dict = self._build_results_dict_trials(rmsds)
        plotter = Plotter(self.prefix, results_dict['results_tuples'])
        plotter.plot_multiple_scatter(data_list=results_dict['results_tuples'],
                                      labels=results_dict['results_labels'],
                                      annotations=self.annot,
                                      colors=results_dict['results_ranks'],
                                      filename='ground_vs_alt1_plddt_ranked_trials',
                                      colorbar_label='Ranking')