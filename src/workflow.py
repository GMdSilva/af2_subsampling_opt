"""
Program main flow
"""

import os
from typing import Dict, Any
import logging

from src.prediction_engine.af2runner import AF2Runner
from src.prediction_engine.msabuilder import MSABuilder
from src.subsampling_optimization.subsamplingoptimizer import SubsamplingOptimizer
from src.mutation_analysis.mutationanalyzer import MutationAnalyzer
from src.utilities.utilities import load_from_pickle
from src.subsampling_optimization.statefinder import StateFinder
from src.mutation_analysis.clustercomparer import ClusterComparer
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
                            '../results',
                            'msas',
                            f"{prefix}_hmmer.a3m")
    predictor = AF2Runner(prefix, msa_path)
    predictor.run_subsampled_af2()


class OptimizerManager:

    def __init__(self, prefix: str):
        self.optimizer = SubsamplingOptimizer(prefix)
        print("Starting Optimizer")

    def get_parameter_set_variations(self):
        print(f"Measuring structural regions of significant variation.")
        significant_ranges = self.optimizer.analyze_predictions('rmsf')
        return self.optimizer.analyze_parameter_set(significant_ranges)

    def get_best_parameter_set(self, analyzed_parameter_set) -> Dict[str, Any]:
        print(f"Ranking parameter sets based on scoring criteria.")
        predictions_path = os.path.join(config.PREDICTION_ROOT,
                                        '../results',
                                        'af2_predictions')
        parameter_set = self.optimizer.make_final_decision(analyzed_parameter_set, predictions_path)
        print(f'chosen parameters: {parameter_set}')
        return parameter_set


def test_mutants(prefix: str) -> None:
    """
    Handles mutation testing for a given prefix and optimized parameters.

    Args:
        prefix (str): The prefix for the mutation analysis.

    Returns:
        None
    """
    print(f"Testing mutants for system {prefix}.")
    analyzer = MutationAnalyzer(prefix)
    results_filename = os.path.join(config.PREDICTION_ROOT,
                                    '../results',
                                    'optimization_results',
                                    f"{prefix}_optimizer_results.pkl")
    optimization_results = load_from_pickle(results_filename)

    filename = f"{prefix}_mut_analysis_results.pkl"
    results_path = os.path.join(config.PREDICTION_ROOT,
                                '../results',
                                'mutant_analysis',
                                filename)

    all_mut_results, mut_data = analyzer.load_or_generate_mut_results(
                                                   results_path,
                                                   optimization_results)
    return analyzer.measure_accuracy(all_mut_results, mut_data)


def get_representative_structures(prefix: str) -> None:
    """
    Handles mutation testing for a given prefix and optimized parameters.

    Args:
        prefix (str): The prefix for the mutation analysis.

    Returns:
        None
    """
    print(f"Getting representative structures for system {prefix}.")
    all_trials = ['16_32', '32_64', '64_128', '128_256', '256_512', '512_1024']
    finder = StateFinder(prefix, all_trials)
    results_filename = os.path.join(config.PREDICTION_ROOT,
                                    '../results',
                                    'optimization_results',
                                    f"{prefix}_optimizer_results.pkl")

    optimization_results = load_from_pickle(results_filename)
    finder.get_refs_and_compare(optimization_results)


def compare_mutation_clusters(prefix):
    print(f"Comparing changes in cluster population for mutations of system {prefix}")
    results_filename = os.path.join(config.PREDICTION_ROOT,
                                    '../results',
                                    'optimization_results',
                                    f"{prefix}_optimizer_results.pkl")
    optimization_results = load_from_pickle(results_filename)
    comparer = ClusterComparer(prefix, optimization_results)
    report = comparer.measure_mutation_effects(measure_accuracy=True)
    print(report)



