""""Functions for editing and saving config file"""

import json
import os
from user_settings.config import PREDICTION_ROOT


def set_config(prefix: str,
               path: str,
               reindex: bool = False,
               first_residue: int = 0,
               build_msa: bool = False,
               run_af2: bool = False,
               optimize_parameters: bool = True,
               jackhmmer_path: str = "",
               colabfoldbatch_path: str = "",
               test_mutants: str = False) -> None:

    config_data = {
        "PREFIX": prefix,
        "PATH": path,
        "REINDEX": reindex,
        "FIRST_RESIDUE": first_residue,
        "BUILD_MSA": build_msa,
        "RUN_AF2": run_af2,
        "OPTIMIZE_PARAMETERS": optimize_parameters,
        "JACKHMMER_PATH": jackhmmer_path,
        "COLABFOLDBATCH_PATH": colabfoldbatch_path,
        "TEST_MUTANTS": test_mutants,
        "CUSTOM_RANGE": [],
    }

    with open("config.json", "w", encoding="utf-8") as config_file:
        json.dump(config_data, config_file, indent=4)


def load_config(filename="config.json"):
    try:
        with open(filename, "r", encoding="utf-8") as config_file:
            config_data = json.load(config_file)
            return config_data
    except FileNotFoundError:
        print(f"{filename} not found.")
        return {}
    except json.JSONDecodeError:
        print(f"Error decoding {filename}. Please ensure it's valid JSON.")
        return {}


def config_template(PREFIX):

    set_config(PREDICTION_ROOT,
               PREFIX,
               os.path.join('results',
                            'predictions',
                            ''),
               reindex=False,
               first_residue=0,
               build_msa=False,
               run_af2=False,
               optimize_parameters=True,
               jackhmmer_path="",
               colabfoldbatch_path="",
               test_mutants=True
               )
