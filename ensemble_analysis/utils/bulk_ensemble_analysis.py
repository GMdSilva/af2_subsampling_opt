"""Measures the structural observables for a range of AF2 prediction sets """

import glob
import os
from typing import Dict, List, Optional

from ensemble_analysis.mdanalysisrunner import MDAnalysisRunner
from utilities.utilities import save_to_pickle


def extract_trial_name(folder: str) -> str:
    """
    Extracts trial name from folder name.

    Parameters:
    - folder: Folder name from which trial name needs to be extracted

    Returns:
    - Extracted trial name
    """
    trial = folder.split('_')
    return f'{trial[-2]}:{trial[-1]}'


def bulk_analysis(prefix: str,
                  path_dirs: str,
                  method: str = 'rmsd',
                  selection: str = 'protein and name CA',
                  ref: Optional[str] = None) -> List[Dict[str, str]]:
    """
    Performs bulk analysis on multiple directories.

    Parameters:
    - prefix: Usually the name of the protein being predicted
    - path_dirs: Directory path where subdirectories are located for analysis
    - method: Observable method for predictions
    - selection: Subset of atoms/features for analysis
    - ref: Path for PDB reference for RMSD calculations
    - save_to_disk: Flag indicating whether to save results to disk

    Returns:
    - List of dictionaries containing analysis results.
    """
    folders = [d for d in glob.glob(f'{path_dirs}/{prefix}') if os.path.isdir(d)]
    all_results = []

    for folder in folders:
        print(f"Analyzing {method} of {prefix} prediction, "
              f"max_seq {folder.split('_')[-2]}, "
              f"with selection: {selection}")
        trial = extract_trial_name(folder)
        analysis_runner = MDAnalysisRunner(folder, selection)
        analysis_results = analysis_runner.analyze_predictions(method=method)

        result_dict = {
            "results": analysis_results['results'],
            'residues': analysis_results['residues'],
            'trial': trial,
            "method": method,
            "selection": selection,
            "reference": ref,
        }
        all_results.append(result_dict)

    filename = os.path.join('results',
                            'mdanalysis_results',
                            f'{prefix}_{method}_results.pkl')
    save_to_pickle(filename, all_results)

    return all_results
