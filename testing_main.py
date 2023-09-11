"""
Program main flow
"""


import os
from typing import Dict, Any

from src.prediction_engine.af2runner import AF2Runner
from src.prediction_engine.msabuilder import MSABuilder
from src.subsampling_optimization.subsamplingoptimizer import SubsamplingOptimizer
from src.mutation_analysis.mutationanalyzer import MutationAnalyzer
from src.utilities.utilities import load_from_pickle
from src.subsampling_optimization.statefinder import StateFinder
from src.mutation_analysis.clustercomparer import ClusterComparer
from user_settings import config


def build_msa(prefix: str) -> None:
    """
    Builds an MSA (Multiple Sequence Alignment) for a given prefix using the MSABuilder.

    Args:
        prefix (str): The prefix for the target protein sequence.

    Returns:
        None
    """
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
    msa_path = os.path.join(config.PREDICTION_ROOT,
                            'results',
                            'msas',
                            f"{prefix}_hmmer.a3m")
    predictor = AF2Runner(prefix, msa_path)
    predictor.run_subsampled_af2()


def optimize_parameters(prefix: str) -> Dict[str, Any]:
    """
    Optimizes parameters for subsampling using the SubsamplingOptimizer.

    Args:
        prefix (str): The prefix for the target protein sequence.

    Returns:
        Dict[str, Any]: A dictionary containing optimized parameters.
    """
    predictions_path = os.path.join(config.PREDICTION_ROOT,
                                    'results',
                                    'af2_predictions')
    method = 'rmsf'
    optimizer = SubsamplingOptimizer(prefix)
    subsampling_results = optimizer.analyze_predictions(method)
    return optimizer.get_optimized_parameters(predictions_path,
                                              subsampling_results)


def test_mutants(prefix: str) -> None:
    """
    Handles mutation testing for a given prefix and optimized parameters.

    Args:
        prefix (str): The prefix for the mutation analysis.

    Returns:
        None
    """
    analyzer = MutationAnalyzer(prefix)
    results_filename = os.path.join(config.PREDICTION_ROOT,
                                    'results',
                                    'optimization_results',
                                    f"{prefix}_optimizer_results.pkl")
    optimization_results = load_from_pickle(results_filename)

    filename = f"{prefix}_mut_analysis_results.pkl"
    results_path = os.path.join(config.PREDICTION_ROOT,
                                'results',
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
    finder = StateFinder(prefix)
    results_filename = os.path.join(config.PREDICTION_ROOT,
                                    'results',
                                    'optimization_results',
                                    f"{prefix}_optimizer_results.pkl")

    optimization_results = load_from_pickle(results_filename)
    finder.get_refs_and_compare(optimization_results)
    comparer = ClusterComparer(prefix, optimization_results)
    report = comparer.measure_mutation_effects(measure_accuracy=True)
    print(report)


def main() -> None:
    """
    Main function to coordinate the pipeline of building MSA,
    running AF2, optimizing parameters, and testing mutants.

    Returns:
        None
    """
    target_prot = 'abl'
    kind = "wt"
    prefix = f"{target_prot}_{kind}"

    if config.BUILD_MSA:
        build_msa(prefix)
    else:
        print("Skipping MSA building, "
              "make sure that path to MSA is accurate in config file")
    if config.RUN_AF2:
        run_af2(prefix)
    else:
        print("Skipping Preliminary AF2 run, "
              "make sure that path to Wild-Type Predictions "
              "is accurate in config file")
    if config.OPTIMIZE_PARAMETERS:
        optimize_parameters(prefix)
    if config.TEST_MUTANTS:
        # muts = load_config('user_settings/mutants.json')
        # for mut in muts:
        #     prefix = mut
        #     build_msa(prefix)
        #     run_af2(prefix)
        accuracy = test_mutants(prefix)
    get_representative_structures(prefix)

if __name__ == "__main__":
    main()
