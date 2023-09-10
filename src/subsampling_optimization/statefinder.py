"""
Defines a class that finds and saves representative structures for each selected peak.
"""

from typing import Dict, Any, Tuple, List
import os
import numpy as np
import glob
import shutil

from src.utilities.plotter import Plotter
from user_settings.new_config import load_config
from user_settings import config
from src.utilities.utilities import load_from_pickle
from src.ensemble_analysis.mdanalysisrunner import MDAnalysisRunner


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
        self.selection: str = 'protein and name CA'
        self.trial: str = None

    def get_variables(self, optimization_results: Dict[str, Any]) -> None:
        """
        Sets trial and range based on given results path and optimization results.

        Args:
            optimization_results (dict): Dictionary containing optimization results.
        """
        self.selection = f"{optimization_results['ranges'][0][0]}_" \
                         f"{optimization_results['ranges'][0][1]}"
        self.trial = optimization_results['parameters'].replace(":", "_")

    def get_reps(self, optimization_results: Dict[str, Any]) -> None:
        """
        Finds the representative structures for each peak.
        """
        self.get_variables(optimization_results)
        peak_calling_path = os.path.join(
            config.PREDICTION_ROOT, 'results', 'peaks',
            f"{self.prefix}_accepted_{self.trial}"
            f"_range_{self.selection}.pkl")
        peak_calling_results = load_from_pickle(peak_calling_path)

        rmsd_path = os.path.join(
            config.PREDICTION_ROOT, 'results', 'optimizer_results',
            f"{self.prefix}_rmsd_{self.trial}_{self.selection}_results.pkl")
        rmsd_results = load_from_pickle(rmsd_path)

        # Finding index for ground
        index_ground = StateFinder._find_closest_index(
            rmsd_results[0]['results'], peak_calling_results['modes'][1])
        prediction_ground = rmsd_results[0]['residues'][index_ground]

        # Finding index for alt
        index_alt = StateFinder._find_closest_index(
            rmsd_results[0]['results'], peak_calling_results['modes'][0])
        prediction_alt1 = rmsd_results[0]['residues'][index_alt]

        return{'ground_ref_index': int(prediction_ground) + 1,
               'alt1_ref_index': int(prediction_alt1) + 1}

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

    def find_pdb_files(self, optimization_results: Dict[str, int]) -> Tuple[List[str], List[str]]:
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

        labeled_files = list(zip(ground_files, ["ground"] * len(ground_files))) + \
                        list(zip(alt1_files, ["alt1"] * len(alt1_files)))

        rmsds = {}
        for file, label in labeled_files:
            self._save_rep_states(file, label)
            rmsd = self.get_rmsd_vs_refs(file, label)
            rmsds[label] = rmsd

        self._plot_results(rmsds)

        return rmsds

    def _construct_path(self, index: int) -> str:
        """
        Constructs a file path based on a given index.

        Args:
            index (int): The index used in the file name.

        Returns:
            str: The constructed file path.
        """
        formatted_index = str(index).zfill(3)
        return os.path.join('results', 'predictions',
                            f"{self.prefix}_{self.trial}",
                            f"{self.prefix}_unrelaxed_rank_{formatted_index}*.pdb")

    def _save_rep_states(self, file, label):
        # Define source and destination paths
        dest_folder = os.path.join('results',
                                   'references',
                                   f"{self.prefix}_"
                                   f"{self.trial}_"
                                   f"{self.selection}_"
                                   f"{label}.pdb")

        # Copy the file
        shutil.copy2(file, dest_folder)

    def get_rmsd_vs_refs(self, file, label):
        ref = file
        prefixes = self._get_prefixes()
        path = os.path.join('results',
                            'predictions')
        results = {}
        rmsds = {}
        trial = self.trial.replace("_", ":")
        selection = f"resid {self.selection.replace('_', ':')} and name CA"


        for prefix in prefixes:
            file_exists = self._check_if_file_exists(prefix, label)
            if file_exists:
                res = file_exists
            else:
                runner = MDAnalysisRunner(path, prefix, selection, ref_name=ref)
                res = runner.process_results(bulk=False,
                                             trial=trial,
                                             method='rmsd_ref',
                                             label=label)
            results[prefix] = res

        return results

    def _check_if_file_exists(self, prefix, label):
        file_exists = os.path.join(config.PREDICTION_ROOT,
                                   'results',
                                   'optimizer_results',
                                   f"{prefix}_rmsd_ref_"
                                   f"{self.trial}_"
                                   f"{self.selection}_"
                                   f"{label}_"
                                   f"results.pkl")
        if os.path.exists(file_exists):
            file = load_from_pickle(file_exists)
            return file
        return False

    @staticmethod
    def _sort_mut_list():
        mut_list = load_config(os.path.join(config.PREDICTION_ROOT,
                                            'user_settings',
                                            'mutants.json'))

        mut_list = dict(sorted(mut_list.items(),
                               key=lambda item: item[1]["rank"]))

        return mut_list

    @staticmethod
    def _get_prefixes():
        mut_list = StateFinder._sort_mut_list()

        prefixes = []
        for mut in mut_list:
            prefixes.append(mut)

        return prefixes

    @staticmethod
    def _get_effects():
        muts_list = StateFinder._sort_mut_list()
        prefixes = StateFinder._get_prefixes()

        effects = []
        for prefix in prefixes:
            effect = muts_list[prefix]['effect']
            effect_label = f"Ground: {str(effect['ground_pop'])} " \
                           f"Alt1: {str(effect['alt1_pop'])}"
            effects.append(effect_label)

        return effects

    @staticmethod
    def _get_mutation_names():
        muts_list = StateFinder._sort_mut_list()
        prefixes = StateFinder._get_prefixes()

        names = []
        for prefix in prefixes:
            name = muts_list[prefix]['label']
            names.append(name)

        return names

    @staticmethod
    def _get_mutation_names():
        muts_list = StateFinder._sort_mut_list()
        prefixes = StateFinder._get_prefixes()

        names = []
        for prefix in prefixes:
            name = muts_list[prefix]['label']
            names.append(name)

        return names

    def _plot_results(self, rmsd_dict):
        labels = StateFinder._get_mutation_names()
        prefixes = StateFinder._get_prefixes()
        effects = StateFinder._get_effects()
        dict = {}
        tupls = []
        for prefix in prefixes:
            dict[f"{prefix}_vs_ground"] = rmsd_dict['ground'][prefix][0]['results']
            dict[f"{prefix}_vs_alt1"] = rmsd_dict['alt1'][prefix][0]['results']
            tupl = (rmsd_dict['ground'][prefix][0]['results'],
                    rmsd_dict['alt1'][prefix][0]['results'])
            tupls.append(tupl)
        plotter = Plotter(self.prefix, tupls)
        plotter.plot_multiple_scatter(tupls, labels, effects)

