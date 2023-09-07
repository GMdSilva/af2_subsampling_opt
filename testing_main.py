"""Scratch file for testing"""

import os
from glob import glob

from user_settings.new_config import load_config
from user_settings import config
from subsampling_optimization.subsamplingoptimizer import SubsamplingOptimizer
from prediction_engine.msabuilder import MSABuilder
from prediction_engine.af2runner import AF2Runner


if __name__ == "__main__":

    target_prot = 'abl'
    kind = "wt"
    prefix = f"{target_prot}_{kind}"

    if config.BUILD_MSA:
        builder = MSABuilder(prefix, sequence="PLACEHLDER")
        builder.build_jackhmmer_msa()

    if config.RUN_AF2:
        msa_path = os.path.join('results', 'msas', '')
        msa_file = glob(os.path.join(msa_path, f"{prefix}_hmmer_*.a3m"))[0]
        predictor = AF2Runner(prefix, msa_file)
        predictor.run_subsampled_af2()

    if config.OPTIMIZE_PARAMETERS:
        method = 'rmsf'
        optimizer = SubsamplingOptimizer(prefix)
        subsampling_results = optimizer.analyze_predictions(config.PATH, method)
        final_ranges = optimizer.get_optimized_parameters(config.PATH, subsampling_results)

    if config.TEST_MUTANTS:
        muts = load_config('user_settings/mutants.json')
        print(muts)

        prefix = f"{target_prot}_{kind}"

