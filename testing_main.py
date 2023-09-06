"""Scratch file for testing"""

import os

from ensemble_analysis.bulk_ensemble_analysis import bulk_analysis
from subsampling_opt.parameter_scoring import final_subsampling_decision
from subsampling_opt.find_ranges import get_final_ranges
from utilities.utilities import load_from_pickle

PATH = os.path.join("examples",
                    "abl",
                    "wild-type")

METHOD = 'rmsf'
PREFIX = 'abl'

RESULTS_PATH = os.path.join("results",
                            "mdanalysis_results",
                            f"{PREFIX}_{METHOD}_results.pkl")

file_exists = os.path.isfile(RESULTS_PATH)
if not file_exists:
    subsampling_results = bulk_analysis(PATH, METHOD)
else:
    subsampling_results = load_from_pickle\
        (f'results/mdanalysis_results/{PREFIX}_{METHOD}_results.pkl')

final_ranges = get_final_ranges(subsampling_results)
subsampling_parameters = final_subsampling_decision(final_ranges, PATH)
