"""
Program main flow
    TODO: REPLACE ABL1 PLACEHOLDER CONSTANTS
        WITH VARIABLES
"""


import os
from glob import glob
from typing import List, Dict, Any

from prediction_engine.af2runner import AF2Runner
from prediction_engine.msabuilder import MSABuilder
from subsampling_optimization.subsamplingoptimizer import SubsamplingOptimizer
from mutation_analysis.mutationanalyzer import MutationAnalyzer
from user_settings import config
from user_settings.new_config import load_config
from utilities.utilities import save_to_pickle, load_from_pickle


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
    msa_path = os.path.join('results', 'msas', '')
    msa_file = glob(os.path.join(msa_path, f"{prefix}_hmmer_*.a3m"))[0]
    predictor = AF2Runner(prefix, msa_file)
    predictor.run_subsampled_af2()


def optimize_parameters(prefix: str) -> Dict[str, Any]:
    """
    Optimizes parameters for subsampling using the SubsamplingOptimizer.

    Args:
        prefix (str): The prefix for the target protein sequence.

    Returns:
        Dict[str, Any]: A dictionary containing optimized parameters.
    """
    method = 'rmsf'
    optimizer = SubsamplingOptimizer(prefix)
    subsampling_results = optimizer.analyze_predictions(method)
    return optimizer.get_optimized_parameters(config.PATH, subsampling_results)


def load_or_generate_mut_results(prefix: str,
                                 results_path: str,
                                 final_ranges: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Loads or generates mutation results for a given prefix and range.

    Args:
        prefix (str): The prefix for the mutation analysis.
        results_path (str): Path where the results are or should be saved.
        final_ranges (Dict[str, Any]): The optimized parameters for subsampling.

    Returns:
        List[Dict[str, Any]]: A list of mutation results.
    """
    old_prefix = prefix
    all_mut_results = []
    file_exists = os.path.isfile(results_path)
    if file_exists:
        return load_from_pickle(results_path)
    muts = load_config('user_settings/mutants.json')
    for mut in muts:
        prefix = mut
        mut_analyzer = MutationAnalyzer(prefix,
                                        final_ranges['ranges'][0]['selection'])
        trial = [[(final_ranges['parameters'])]]
        mut_results = mut_analyzer.compare_conditions(trials=trial,
                                                      old_prefix=old_prefix)
        all_mut_results.append(mut_results)
    save_to_pickle(results_path, all_mut_results)
    return all_mut_results


def plot_mut_results(prefix: str,
                     all_mut_results: List[Dict[str, Any]],
                     final_ranges: Dict[str, Any]) -> None:
    """
    Plots mutation results for a given prefix and range.

    Args:
        prefix (str): The prefix for the mutation analysis.
        all_mut_results (List[Dict[str, Any]]): A list of mutation results.
        final_ranges (Dict[str, Any]): The optimized parameters for subsampling.

    Returns:
        None
    """
    mut_analyzer = MutationAnalyzer(prefix, final_ranges['ranges'][0]['selection'])
    labels = ['ground_pop_diff', 'alt1_pop_diff', 'ground_pop_test', 'alt1_pop_test']
    for label in labels:
        mut_analyzer.plot_results(all_mut_results, label)


def test_mutants(prefix: str, final_ranges: Dict[str, Any]) -> None:
    """
    Handles mutation testing for a given prefix and optimized parameters.

    Args:
        prefix (str): The prefix for the mutation analysis.
        final_ranges (Dict[str, Any]): The optimized parameters for subsampling.

    Returns:
        None
    """
    filename = "abl_results"
    results_path = os.path.join('results', 'mutation_analysis', filename)
    all_mut_results = load_or_generate_mut_results(prefix, results_path, final_ranges)
    plot_mut_results(prefix, all_mut_results, final_ranges)


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
    if config.RUN_AF2:
        run_af2(prefix)
    if config.OPTIMIZE_PARAMETERS:
        final_ranges = optimize_parameters(prefix)
    if config.TEST_MUTANTS:
        # muts = load_config('user_settings/mutants.json')
        # for mut in muts:
        #     prefix = mut
        #     build_msa(prefix)
        #     run_af2(prefix)
        test_mutants(prefix, final_ranges)


if __name__ == "__main__":
    main()
