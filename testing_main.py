"""Scratch file for testing"""

import os

from user_settings.new_config import set_config
from analyze_predictions import analyze_predictions, get_optimized_parameters



PATH = os.path.join("examples",
                    "abl",
                    "wild-type")

METHOD = 'rmsf'
PREFIX = 'abl'

set_config(PREFIX, PATH)

subsampling_results = analyze_predictions(PATH, METHOD)
final_ranges = get_optimized_parameters(PATH, subsampling_results)
