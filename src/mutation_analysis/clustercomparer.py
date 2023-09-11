"""
Module for finding and saving representative structures for each selected peak.
"""

from typing import List, Dict, Union, Tuple, Optional
import glob
import os
import warnings

import json
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from kneed import KneeLocator

from user_settings.config import SYSTEM_NAME, IS_JUPYTER
from user_settings.new_config import load_config
from user_settings.config import PREDICTION_ROOT
from src.utilities.plotter import Plotter
from src.mutation_analysis.mutantstatefinder import MutantStateFinder


warnings.filterwarnings("ignore",
                        category=UserWarning,
                        module='sklearn',
                        message='.*memory leak on Windows with MKL.*')


class ClusterComparer(MutantStateFinder):
    """
    Analyzes mutations and measures accuracy at predicting the effects of mutations in
    the ground or alternative states.
    """

    def __init__(self, prefix: str, optimization_results: dict) -> None:
        """
        Initialize the ClusterComparer instance.

        Args:
        - prefix (str): Usually the name of the protein being studied.
        - optimization_results (dict): Results obtained after optimization.
        """
        self.prefix = prefix
        self.selection: str = 'protein and name CA'
        self.wildtype_results: dict = None
        self.rmsd_dict = self.get_refs_and_compare_muts(optimization_results)

    def get_wildtype_results(self) -> np.ndarray:
        """
        Get wildtype results as a 2D array.

        Returns:
        - np.ndarray: A 2D numpy array containing wildtype results.
        """
        wildtype_tuple = (
            self.rmsd_dict['ground'][SYSTEM_NAME][0]['results'],
            self.rmsd_dict['alt1'][SYSTEM_NAME][0]['results']
        )
        return np.column_stack(wildtype_tuple)

    @staticmethod
    def find_optimal_clusters(wildtype_results: np.ndarray,
                              plot_results: bool = True,
                              save_results: bool = True) -> int:
        """
        Find optimal number of clusters using the elbow method.

        Parameters:
        - wildtype_results (np.ndarray): A 2D numpy array containing wildtype results.
        - plot_results (bool, optional): If True, plots the elbow curve. Defaults to True.
        - save_results (bool, optional): If True, saves the elbow curve plot. Defaults to True.

        Returns:
        - int: Optimal number of clusters determined by the elbow method.

        Note:
        The function uses the distortion as a metric to determine the optimal number of clusters.
        It also uses the "kneed" library to automatically detect the elbow point. If the elbow
        point isn't detected, a default value of 3 clusters is returned.
        """
        distortions = []
        k_tests = range(1, 25)

        for k_test in k_tests:
            if os.name == 'nt':
                # Set the environment variable
                os.environ['OMP_NUM_THREADS'] = '1'
            kmean_model = KMeans(n_clusters=k_test, n_init=10)
            kmean_model.fit(wildtype_results)
            distortions.append(kmean_model.inertia_)

        if plot_results:
            plt.figure(figsize=(8, 6))
            plt.plot(k_tests, distortions, 'bo-')
            plt.xlabel('Number of Clusters')
            plt.ylabel('Distortion')
            plt.title(f'Elbow Method showing optimal k for {SYSTEM_NAME}')
            if save_results:
                plt.savefig(os.path.join(PREDICTION_ROOT,
                                         'results',
                                         'plots',
                                         f"{SYSTEM_NAME}_kmeans_elbow_test.png"))
            if IS_JUPYTER:
                plt.show()

        # Using kneed library to find the elbow point
        knee = KneeLocator(k_tests, distortions, curve='convex', direction='decreasing')
        k_optimal = knee.elbow

        if k_optimal is None:  # Just in case the kneed method can't find the elbow
            k_optimal = 3

        return k_optimal

    @staticmethod
    def evaluate_clusters(cluster_labels: np.ndarray) -> np.ndarray:
        """
        Evaluate and return counts of data points in each cluster.

        Given the labels of data points indicating their respective clusters,
        this function computes and returns the number of data points in each cluster.

        Parameters:
        - cluster_labels (np.ndarray): An array containing cluster assignments for each data point.

        Returns:
        - np.ndarray: An array containing counts of data points for each unique cluster.

        Example:
        If cluster_labels = [0, 1, 0, 2, 1], the function will return [2, 2, 1],
        meaning there are 2 data points in cluster 0, 2 in cluster 1, and 1 in cluster 2.
        """
        _, counts = np.unique(cluster_labels, return_counts=True)
        return counts

    def measure_mutation_effects(self, measure_accuracy: bool = True) -> None:
        """
        Measure mutation effects on cluster populations and plot the results.

        Parameters:
        - measure_accuracy (bool): Whether to measure accuracy or not. Default is True.

        Returns:
        - None
        """
        # Build wildtype model and get RMSD measurements
        kmeans, clusters_wildtype = self.build_wt_model()
        all_rmsd_measurements = self._build_results_dict(self.rmsd_dict)

        # Get clustering results
        clustering_results = self._get_clustering_results(all_rmsd_measurements,
                                                          clusters_wildtype,
                                                          kmeans)

        # Plot cluster results
        self._plot_cluster_data(all_rmsd_measurements, clustering_results)

        # Generate reports, with optional accuracy measurement
        report = self._generate_report(measure_accuracy,
                                       all_rmsd_measurements,
                                       clustering_results)
        return report

    def build_wt_model(self) -> Tuple[KMeans, np.ndarray]:
        """
        Build a wildtype model and determine the optimal clusters for the wildtype results.

        This function retrieves the wildtype results and uses the elbow method to
        determine the optimal number of clusters. It then builds a KMeans clustering model
        using this optimal number and returns both the trained model and the resulting clusters.

        Returns:
        - kmeans (KMeans): The trained KMeans clustering model.
        - clusters_wildtype (np.ndarray): An array of cluster assignments for each data point.

        Note:
        The method uses the wildtype results and aims to find an optimal cluster number
        using the elbow method. It then trains the KMeans clustering model using this optimal
        number and assigns each data point to a specific cluster.
        """
        wildtype_results = self.get_wildtype_results()
        k_optimal = self.find_optimal_clusters(wildtype_results)
        kmeans = KMeans(n_clusters=k_optimal)
        clusters_wildtype = kmeans.fit_predict(wildtype_results)

        return kmeans, clusters_wildtype

    def _plot_cluster_data(self, all_rmsd_measurements: dict, clustering_results: dict) -> None:
        """
        Plot cluster results using the provided measurement and clustering data.

        Parameters:
        - all_rmsd_measurements (dict): Dictionary containing all RMSD measurement data.
        - clustering_results (dict): Dictionary containing clustering data
                like cluster colors and annotations.

        Returns:
        - None. The function will internally call another function to handle plotting.
        """

        self.plot_cluster_results(all_rmsd_measurements,
                                  clustering_results['cluster_colors'],
                                  clustering_results['cluster_annotations'])

    def _generate_report(self,
                         measure_accuracy: bool,
                         all_rmsd_measurements: dict,
                         clustering_results: dict) -> None:
        """
        Generate a reports based on the clustering results and optionally measure accuracy.

        Parameters:
        - measure_accuracy (bool): A flag to determine if accuracy should be measured or not.
        - all_rmsd_measurements (dict): Dictionary containing all RMSD measurement data.
        - clustering_results (dict): Dictionary containing clustering results and evaluations.

        Returns:
        - None. The function will internally call another function to handle reports generation.
        """

        if measure_accuracy:
            accuracy = self.measure_accuracy(all_rmsd_measurements['results_labels'],
                                             all_rmsd_measurements['results_effects'],
                                             clustering_results['results'])

            report = self.build_report(clustering_results['wt_pct'],
                                       clustering_results['results'],
                                       accuracy_report=accuracy)
        else:
            report = self.build_report(clustering_results['wt_pct'],
                                       clustering_results['results'])
        return report

    def _get_clustering_results(self,
                                all_rmsd_measurements: Dict[str, List],
                                clusters_wildtype: List,
                                kmeans: KMeans) -> Dict[str, Union[np.ndarray, List]]:
        """
        Helper method to get clustering results.

        Parameters:
        - all_rmsd_measurements (dict): Dictionary of RMSD measurements.
        - clusters_wildtype (list): List of wildtype clusters.
        - kmeans (KMeans): KMeans model instance.

        Returns:
        - Dictionary containing clustering results, cluster colors, and cluster annotations.
        """
        cluster_annotations = []
        cluster_colors = []
        clustering_results = {}

        wt_evaluation = self.evaluate_clusters(clusters_wildtype)
        sorted_indices = np.argsort(wt_evaluation)
        wt_evaluation = wt_evaluation[sorted_indices][::-1]
        _, wt_pct = self._format_cluster_populations(wt_evaluation)
        wt_cluster_centroids = kmeans.cluster_centers_
        wt_centroid_indices = self._get_centroid_indices(wt_cluster_centroids[sorted_indices])
        centroid_labels = []
        for i in range(0, len(wt_centroid_indices)):
            label = f"s{i+1}"
            centroid_labels.append(label)
        self._save_centroids(wt_centroid_indices, centroid_labels)

        for result_tuple, label in zip(all_rmsd_measurements['results_tuples'],
                                       all_rmsd_measurements['results_labels']):
            clusters_distribution = kmeans.predict(np.column_stack(result_tuple))

            mut_evaluation = self._evaluate_mutations(clusters_distribution,
                                                      sorted_indices,
                                                      wt_evaluation)

            clustering_results[label] = [mut_evaluation['clusters_diff'],
                                         mut_evaluation['all_clusters']]
            cluster_colors.append(clusters_distribution)

            formatted, _ = self._format_cluster_populations(mut_evaluation['all_clusters'])
            cluster_annotations.append(formatted)

        return {
            'results': clustering_results,
            'cluster_colors': cluster_colors,
            'cluster_annotations': cluster_annotations,
            'wt_evaluation': wt_evaluation,
            'wt_pct': wt_pct
        }

    def _evaluate_mutations(self, clusters_distribution1: np.ndarray,
                            sorted_indices: np.ndarray,
                            wt_evaluation: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Helper method to evaluate mutations.

        Parameters:
        - clusters_distribution1 (np.ndarray): The cluster distribution of the mutations.
        - sorted_indices (np.ndarray): Indices for sorting clusters.
        - wt_evaluation (np.ndarray): Wildtype evaluation array.

        Returns:
        - Dictionary containing clusters difference and all clusters.
        """
        mut_evaluation = self.evaluate_clusters(clusters_distribution1)
        mut_evaluation = mut_evaluation[sorted_indices][::-1]
        clusters_diff = mut_evaluation - wt_evaluation
        return {
            'clusters_diff': clusters_diff,
            'all_clusters': mut_evaluation
        }

    @staticmethod
    def load_and_sort_mut_data() -> Tuple[Dict[str, dict], List[str]]:
        """
        Load and sort mutation data from a predefined configuration.

        Returns:
            Tuple containing sorted mutation data and mutation data names.
        """
        mut_data = load_config('user_settings/mutants.json')
        sorted_mut_data = dict(sorted(mut_data.items(), key=lambda item: item[1]["rank"]))
        mut_data_names = [data['label'] for data in sorted_mut_data.values()]
        return sorted_mut_data, mut_data_names

    def compute_results(self, clustering_results: Dict[str, Tuple[List[float], List[float]]]) \
            -> Tuple:
        """
        Compute total and differential results based on clustering results.

        Args:
            clustering_results: Dictionary containing clustering results.

        Returns:
            A tuple containing dictionaries of total and differential results,
            lists of states and labels, and counts of ground up and ground down values.
        """
        total_results = {}
        diff_results = {}
        states = set()
        labels_up = []
        labels_down = []
        count_ground_up = 0
        count_ground_down = 0

        # Assuming every label has the same number of states
        total_measurements = sum(next(iter(clustering_results.values()))[1])

        for label, (diff_vals, total_vals) in clustering_results.items():
            total_results[label] = self.calculate_population(total_vals,
                                                             total_measurements)
            diff_res, ground_up, ground_down =\
                self.calculate_population_difference(diff_vals, total_measurements)

            if ground_up:
                count_ground_up += 1
                labels_up.append(label)
            if ground_down:
                count_ground_down += 1
                labels_down.append(label)

            diff_results[label] = diff_res
            states.update(diff_res.keys())

        return total_results,\
               diff_results,\
               list(states),\
               labels_up,\
               labels_down,\
               count_ground_up,\
               count_ground_down

    @staticmethod
    def calculate_population(values: List[float], total_measurements: float) -> Dict[str, float]:
        """
        Calculate population percentages from the given values.

        Args:
            values: List of population values.
            total_measurements: Total measurement value.

        Returns:
            Dictionary with state labels and their respective population percentages.
        """
        return {f"S{s}" if s else "Ground":
                    (val / total_measurements) * 100 for s, val in enumerate(values)}

    @staticmethod
    def calculate_population_difference(values: List[float],
                                        total_measurements: float) -> Tuple:
        """
        Calculate population differences from the given values.

        Args:
            values: List of population values.
            total_measurements: Total measurement value.

        Returns:
            Tuple containing dictionary with state labels and their respective population differences, and
            two booleans indicating ground up and ground down values.
        """
        ground_up, ground_down = False, False
        results = {}
        for s, val in enumerate(values):
            state = f"S{s}" if s else "Ground"
            results[state] = (val / total_measurements) * 100

            if s == 0:
                ground_up = val > 0
                ground_down = val < 0

        return results, ground_up, ground_down

    def build_report(self, wt_evaluation: Union[List[str], Tuple[str, ...]],
                     clustering_results: Dict[str, Tuple[List[float], List[float]]],
                     accuracy_report: Optional[Dict[str, Union[float, List[str]]]] = None) -> None:
        """
        Build and print a reports based on clustering results, and optionally an accuracy reports.

        Args:
            wt_evaluation: Evaluations for the wild-type data.
            clustering_results: Dictionary containing clustering results and evaluations.
            accuracy_report: Optional dictionary containing accuracy measurements.

        Returns:
            None. The reports is printed to the console.
        """
        sorted_mut_data,\
        mut_data_names = self.load_and_sort_mut_data()
        total_results,\
        diff_results,\
        states,\
        labels_up,\
        labels_down,\
        count_ground_up,\
        count_ground_down = self.compute_results(clustering_results)

        plotter = Plotter(SYSTEM_NAME, diff_results)
        for state in states:
            plotter.plot_bar_for_state(state, sorted_mut_data)

        report_dict = {
            'System Name': SYSTEM_NAME,
            'Optimized max_seq': self.trial.split('_')[0],
            'Optimized extra_seq': self.trial.split('_')[1],
            'C-Alpha Selection': f"{self.selection.split('_')[0]}:"
                                 f"{self.selection.split('_')[1]}",
            'Detected WT States': len(wt_evaluation),
            'WT Ground Population': wt_evaluation[0].split(':')[1],
            'Alternative State Populations': wt_evaluation[1:],
            '# Tested Mutations': len(sorted_mut_data)-1,
            'Tested Variants': mut_data_names,
            '# Ground State Stabilizing Variants': count_ground_up,
            '# Ground State Destabilizing Variants': count_ground_down,
            'Ground State Stabilizing Variants': labels_up,
            'Ground State Destabilizing Variants': labels_down,
            'Relative State Populations': total_results,
            'Relative State Population Differences': diff_results
        }

        # Extend the reports if accuracy reports is given
        if accuracy_report:
            report_dict.update({
                'Accuracy %': accuracy_report['accuracy_%'],
                'Correct Predictions': accuracy_report['correct_predictions'],
                'Incorrect Predictions': accuracy_report['failed_predictions']
            })
        report_path = os.path.join('results',
                                   'reports',
                                   f"{SYSTEM_NAME}_results_report.json")
        # Writing JSON data to a file
        with open(report_path, 'w') as file:
            json.dump(report_dict, file)
        return report_dict

    @staticmethod
    def measure_accuracy(labels: List[str],
                         expected_effects: List[str],
                         clustering_results: Dict[str, Union[int, List[int]]]) -> \
            Dict[str, Union[float, List[List[str]]]]:
        """
        Measures the accuracy of clustering results against the expected effects.

        Parameters:
        - labels (List[str]): List of labels for the measurements.
        - expected_effects (List[str]): List of expected effects, each containing a string
          representation like 'Ground: +'.
        - clustering_results (Dict[str, Union[int, List[int]]]): Dictionary containing clustering
          results. Each key corresponds to a label and the value contains the result.

        Returns:
        - Dict[str, Union[float, List[List[str]]]]: A dictionary containing the accuracy percentage,
          correct predictions, and failed predictions.
        """

        total_measurements = len(expected_effects)
        ground_truths = [1 if 'Ground: +' in effect else (-1 if 'Ground: -' in effect else 'r')
                         for effect in expected_effects]

        results = [clustering_results[label][0][0] for label in labels]

        correct_predictions, failed_predictions = [], []
        accuracy_count = sum(1 for ground_truth, result, label in zip(ground_truths,
                                                                      results,
                                                                      labels)
                             if ground_truth == int(np.sign(result)))

        for ground_truth, result, label in zip(ground_truths, results, labels):
            if ground_truth != 'r':
                predicted_effect = '+' if result > 0 else '-'
                actual_effect = '+' if ground_truth > 0 else '-'
                if ground_truth == int(np.sign(result)):
                    correct_predictions.append([label, actual_effect, predicted_effect])
                else:
                    failed_predictions.append([label, actual_effect, predicted_effect])

        accuracy_percentage = (accuracy_count / (total_measurements - 1)) * 100

        accuracy_report = {
            'accuracy_%': accuracy_percentage,
            'correct_predictions': correct_predictions,
            'failed_predictions': failed_predictions,
        }
        return accuracy_report

    @staticmethod
    def _format_cluster_populations(mut_evaluations: np.ndarray) ->\
            Tuple[str, List[str]]:
        """
        Format cluster populations into percentage representation.

        Parameters:
        - mut_evaluations (np.ndarray): Array containing the evaluation values for mutations.

        Returns:
        - Tuple[str, List[str]]: The joined formatted populations
            and a list of formatted population values.
        """
        values_pct = []
        total_samples = mut_evaluations.sum()
        for i, value in enumerate(mut_evaluations):
            value_pct = (value / total_samples) * 100
            value_pct = round(value_pct * 10) / 10
            values_pct.append(str(f"S{i + 1}: {value_pct}%"))
        values_pct_joined = ' '.join(values_pct)
        return values_pct_joined, values_pct

    @staticmethod
    def plot_cluster_results(all_rmsd_measurements: dict,
                             cluster_colors: List[str],
                             annotations: List[str]) -> None:
        """
        Plot cluster results using the provided data.

        Parameters:
        - all_rmsd_measurements (dict): Dictionary containing measurements related to RMSD.
        - cluster_colors (List[str]): List of colors representing different clusters.
        - annotations (List[str]): List of annotations for the clusters.

        Returns:
        - None: The function performs plotting and does not return any value.
        """
        plotter = Plotter(SYSTEM_NAME, all_rmsd_measurements)
        plotter.plot_multiple_scatter(data_list=all_rmsd_measurements['results_tuples'],
                                      labels=all_rmsd_measurements['results_labels'],
                                      annotations=annotations,
                                      colors=cluster_colors,
                                      fontsize=7,
                                      filename='mutations_clustered_wt')

    def _get_centroid_indices(self, cluster_centroids: np.ndarray) -> np.ndarray:
        """
        Get indices of data points closest to each centroid.

        Parameters:
        - cluster_centroids (np.ndarray): Array containing cluster centroids.

        Returns:
        - np.ndarray: Array of indices representing data points closest to each centroid.
        """
        wt_measurements = self.get_wildtype_results()
        distances = np.linalg.norm(wt_measurements[:, np.newaxis] - cluster_centroids, axis=2)
        return np.argmin(distances, axis=0)

    def _save_centroids(self, centroid_indexes: np.ndarray, states: List[str]) -> None:
        """
        Save representative states for given centroid indices.

        Parameters:
        - centroid_indexes (np.ndarray): Indices of centroids.
        - states (List[str]): List of states corresponding to centroids.

        Returns:
        - None: The function saves the states and does not return any value.
        """
        for index, state in zip(centroid_indexes, states):
            if state == 'S1':
                state = 'ground_2d'
            path = self._construct_path(index)
            path = glob.glob(path)
            self._save_rep_states(path[0], state)
