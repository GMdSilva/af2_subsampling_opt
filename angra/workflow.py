"""
Program main flow
"""

import os
from typing import Dict, Any
import logging

from angra.prediction_engine.af2runner import AF2Runner
from angra.prediction_engine.msabuilder import MSABuilder
from angra.subsampling_optimization.subsamplingoptimizer import SubsamplingOptimizer
from angra.mutation_analysis.mutationanalyzer import MutationAnalyzer
from angra.utilities.utilities import load_from_pickle
from angra.subsampling_optimization.statefinder import StateFinder
from angra.mutation_analysis.clustercomparer import ClusterComparer
from user_settings import config

def setup_logging():
    # Set up basic logging configuration
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')


def build_msa(prefix: str) -> None:
    """
    Builds an MSA (Multiple Sequence Alignment) for a given prefix using the MSABuilder.

    Args:
        prefix (str): The prefix for the target protein sequence.

    Returns:
        None
    """
    print("If jackhmmer is in path, building MSA.")
    builder = MSABuilder(prefix, sequence="PLACEHLDER")
    builder.build_jackhmmer_msa()


def run_af2(prefix: str) -> None:
    """
    Runs AF2 (AlphaFold 2) for a given prefix.

    Args:
        prefix (str): The prefix for the target protein sequence.

    Returns:
        None
    """
    print("If AlphaFold2 is in path, running predictions.")
    msa_path = os.path.join(config.PREDICTION_ROOT,
                            'results',
                            'msas',
                            f"{prefix}_hmmer.a3m")
    predictor = AF2Runner(prefix, msa_path)
    predictor.run_subsampled_af2()


class OptimizerManager:

    def __init__(self, prefix: str):
        self.optimizer = SubsamplingOptimizer(prefix)

    def get_variation_regions(self) -> Dict[str, Any]:
        self.optimizer.analyze_predictions('rmsf')

    def plot_variation_regions(self):
        self.optimizer.plot_rmsf_results()

    def get_parameter_set_variations(self):
        self.optimizer.analyze_parameter_set()

    def plot_parameter_set(self):
        self.optimizer.plot_rmsd_ranges()

    def get_best_parameter_set(self) -> Dict[str, Any]:
        predictions_path = os.path.join(config.PREDICTION_ROOT,
                                        'results',
                                        'af2_predictions')
        best_set = self.optimizer.make_final_decision(predictions_path)
        return best_set


class MutantTester:
    def __init__(self, prefix: str):
        """
        Initializes the MutantTester object with the provided prefix.

        Args:
            prefix (str): The prefix for the mutation analysis.
        """
        self.prefix = prefix
        self.results_filename = os.path.join(config.PREDICTION_ROOT,
                                             'results',
                                             'optimization_results',
                                             f"{self.prefix}_optimizer_results.pkl")
        self.optimization_results = None
        self.analyzer = MutationAnalyzer(prefix)
        self.mut_results_filename = f"{self.prefix}_mut_analysis_results.pkl"
        self.results_path = os.path.join(config.PREDICTION_ROOT,
                                         'results',
                                         'mutant_analysis',
                                         self.mut_results_filename)

    def load_optimizer_results(self):
        """Loads optimization results."""
        self.optimization_results = load_from_pickle(self.results_filename)

    def load_or_generate_mut_results(self):
        """
        Load or generate mutation results.

        Returns:
            all_mut_results, mut_data
        """
        return self.analyzer.load_or_generate_mut_results(self.results_path,
                                                          self.optimization_results)

    def test_mutants(self) -> None:
        """
        Handles mutation testing for the initialized prefix and optimized parameters.
        """
        print(f"Testing mutants for system {self.prefix}.")
        self.load_optimizer_results()
        self.load_or_generate_mut_results()
        return self.analyzer.measure_accuracy()

    def plot_results(self):
        self.analyzer.plot_mut_results()


def get_representative_structures(prefix: str, all_trials: list) -> None:
    """
    Handles mutation testing for a given prefix and optimized parameters.

    Args:
        prefix (str): The prefix for the mutation analysis.

    Returns:
        None
    """
    print(f"Getting representative structures for system {prefix}.")
    finder = StateFinder(prefix, all_trials)
    results_filename = os.path.join(config.PREDICTION_ROOT,
                                    'results',
                                    'optimization_results',
                                    f"{prefix}_optimizer_results.pkl")

    optimization_results = load_from_pickle(results_filename)
    finder.get_refs_and_compare(optimization_results)


class MutationClusterManager:
    def __init__(self, prefix):
        self.prefix = prefix
        self.results_filename = os.path.join(config.PREDICTION_ROOT,
                                             'results',
                                             'optimization_results',
                                             f"{self.prefix}_optimizer_results.pkl")
        self.optimization_results = None
        self.comparer = None
        self.kmeans = None
        self.clusters_wildtype = None
        self.report = None


    def build_wt_model(self):
        self.optimization_results = load_from_pickle(self.results_filename)
        self.comparer = ClusterComparer(self.prefix, self.optimization_results)
        self.kmeans, self.clusters_wildtype = self.comparer.build_wt_model()

    def measure_effects(self):
        self.report = self.comparer.measure_mutation_effects(self.kmeans,
                                                             self.clusters_wildtype,
                                                             measure_accuracy=True)

    def plot_cluster_results(self):
        self.comparer.f_plot_cluster_results()

    def plot_clustering_diffs(self):
        self.comparer.f_plot_clustering_diffs()

    def plot_mutant_state_finder_results(self):
        self.comparer.f_plot_mutant_state_finder_results()
