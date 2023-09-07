"""Scratch file for testing"""

import os
from glob import glob

from prediction_engine.af2runner import AF2Runner
from prediction_engine.msabuilder import MSABuilder
from subsampling_optimization.subsamplingoptimizer import SubsamplingOptimizer
from user_settings import config
from user_settings.new_config import load_config


if __name__ == "__main__":

    TARGET_PROT = 'abl'
    KIND = "wt"
    prefix = f"{TARGET_PROT}_{KIND}"

    if config.BUILD_MSA:
        builder = MSABuilder(prefix, sequence="PLACEHLDER")
        builder.build_jackhmmer_msa()

    if config.RUN_AF2:
        msa_path = os.path.join('results', 'msas', '')
        msa_file = glob(os.path.join(msa_path, f"{prefix}_hmmer_*.a3m"))[0]
        predictor = AF2Runner(prefix, msa_file)
        predictor.run_subsampled_af2()

    if config.OPTIMIZE_PARAMETERS:
        METHOD = 'rmsf'
        optimizer = SubsamplingOptimizer(prefix)
        subsampling_results = optimizer.analyze_predictions(METHOD)
        final_ranges = optimizer.get_optimized_parameters(config.PATH, subsampling_results)

    if config.TEST_MUTANTS:
        all_mut_results = []
        old_prefix = prefix
        muts = load_config('user_settings/mutants.json')
        mut_path = os.path.join('results',
                                'predictions',
                                '')
        selection = final_ranges['ranges'][0]['selection']
        trial = [[(final_ranges['parameters'])]]
        for mut in muts:
            prefix = mut
            optimizer = SubsamplingOptimizer(prefix, use_custom_range=True, selection=selection)
            mut_results = optimizer.compare_conditions(trials=trial)
            all_mut_results.append(mut_results)
