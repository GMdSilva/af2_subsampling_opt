"""Measures the structural observables for a range of AF2 prediction sets """

from typing import List, Dict, Optional

import glob
import os

from utilities.utilities import save_to_pickle
from ensemble_analysis.ensemble_analysis import analyze_predictions
from user_settings.config import PREFIX


def extract_trial_name(folder: str) -> str:
    """
    Extracts trial name from folder name.

    Parameters:
    - folder: Folder name from which trial name needs to be extracted

    Returns:
    - Extracted trial name
    """
    trial = folder.split('_')
    return f'{trial[-2]}:{trial[-1][:-1]}'


def bulk_analysis(path_dirs: str,
                  method: str = 'rmsd',
                  selection: str = 'protein and name CA',
                  ref: Optional[str] = None,
                  save_to_disk: bool = True) -> List[Dict[str, str]]:
    """
    Performs bulk analysis on multiple directories.

    Parameters:
    - path_dirs: Directory path where subdirectories are located for analysis
    - method: Observable method for predictions
    - selection: Subset of atoms/features for analysis
    - ref: Path for PDB reference for RMSD calculations
    - save_to_disk: Flag indicating whether to save results to disk

    Returns:
    - List of dictionaries containing analysis results.
    """
    folders = [d for d in glob.glob(f'{path_dirs}/*') if os.path.isdir(d)]
    all_results = []

    for folder in folders:
        print(f"Analyzing {method} of {PREFIX} prediction, "
              f"max_seq {folder.split('_')[-2]}, "
              f"with selection: {selection}")
        trial = extract_trial_name(folder)
        res = analyze_predictions(folder, method, selection)

        result_dict = {
            "results": res['results'],
            'residues': res['residues'],
            'trial': trial,
            "method": method,
            "selection": selection,
            "reference": ref
        }
        all_results.append(result_dict)

    if save_to_disk:
        filename = os.path.join('results',
                                'mdanalysis_results',
                                f'{PREFIX}_{method}_results.pkl')
        save_to_pickle(filename, all_results)

    return all_results
