"""Scratch file for testing"""

import os
from glob import glob

from subsampling_opt.analyze_predictions import analyze_predictions, get_optimized_parameters
from prediction_engine.build_msa import build_jackhmmer_msa
from prediction_engine.make_predictions import run_subsampled_af2


if __name__ == "__main__":

    from user_settings import config
    if config.BUILD_MSA:
        build_jackhmmer_msa('PLACEHLDER')

    if config.RUN_AF2:
        msa_path = os.path.join('results', 'msas', '')
        msa_file = glob(os.path.join(msa_path, f"{config.PREFIX}_hmmer_*.a3m"))[0]
        run_subsampled_af2(msa_file)

    if config.OPTIMIZE_PARAMETERS:
        method = 'rmsf'
        subsampling_results = analyze_predictions(config.PATH, method)
        final_ranges = get_optimized_parameters(config.PATH, subsampling_results)
