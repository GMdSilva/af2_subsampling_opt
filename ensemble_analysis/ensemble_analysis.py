"""Runs MDAnalysis to measure structural observables in AF2 prediction ensembles"""

from typing import Dict

import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import numpy as np

from user_settings.config import REFERENCE_PATH
from utilities.utilities import load_trajectory


def load_rmsd(prediction_results: mda.Universe,
              selection: str = 'protein and name CA') -> Dict[str, np.ndarray]:
    """
    Calculates the Root Mean Square Deviation (RMSD) for given predictions.
    Uses the #1 ranked prediction (by pLDDT) as reference.

    Parameters:
    - prediction_results: Loaded trajectory data
    - selection: Subset of atoms/features to run the analysis on

    Returns:
    - Dictionary with RMSD analysis results.
    """
    pred_rmsd = rms.RMSD(prediction_results,
                         select='protein and name CA',
                         groupselections=[selection])
    pred_rmsd.run()
    return {'residues': pred_rmsd.results.rmsd[:, 1], 'results': pred_rmsd.results.rmsd[:, 3]}


def load_rmsd_ref(prediction_results: mda.Universe,
                  selection: str = 'protein and name CA') -> Dict[str, np.ndarray]:
    """
    Calculates the RMSD using a reference for given predictions.

    Parameters:
    - prediction_results: Loaded trajectory data
    - selection: Subset of atoms/features to run the analysis on

    Returns:
    - Dictionary with RMSD analysis results using a reference.
    """
    reference_path = REFERENCE_PATH
    ref = mda.Universe(reference_path)
    ref = ref.select_atoms(selection)
    pred_rmsd = rms.RMSD(prediction_results,
                         ref,
                         select='protein and name CA',
                         groupselections=[selection])
    pred_rmsd.run()
    return {'residues': pred_rmsd.results.rmsd[:, 1], 'results': pred_rmsd.results.rmsd[:, 3]}


def load_rmsf(prediction_results: mda.Universe,
              selection: str = 'protein and name CA') -> Dict[str, np.ndarray]:
    """
    Calculates the Root Mean Square Fluctuation (RMSF) for given predictions.

    Parameters:
    - prediction_results: Loaded trajectory data
    - selection: Subset of atoms/features to run the analysis on

    Returns:
    - Dictionary with RMSF analysis results.
    """
    average = align.AverageStructure(prediction_results,
                                     prediction_results,
                                     select=selection,
                                     ref_frame=0).run()
    ref = average.results.universe
    align.AlignTraj(prediction_results, ref, select=selection, in_memory=True).run()
    c_alphas = prediction_results.select_atoms(selection)
    pred_rmsf = rms.RMSF(c_alphas).run()
    return {'residues': c_alphas.resids, 'results': pred_rmsf.results.rmsf}


def analyze_predictions(path: str,
                        method: str = 'rmsd',
                        selection: str = 'protein and name CA') -> Dict[str, np.ndarray]:
    """
    Analyzes a trajectory based on the given configuration.

    Parameters:
    - path: Output directory of AlphaFold2 run
    - method: observable to measure in predictions (e.g. RMSD, distance between X and Y, etc.)
    - selection: subset of atoms/features to run analysis on

    Returns:
    - Dictionary containing analysis results.
    """
    prediction_results = load_trajectory(path)

    analysis_methods = {
        'rmsd': load_rmsd,
        'rmsd_ref': load_rmsd_ref,
        'rmsf': load_rmsf
    }

    if method in analysis_methods:
        return analysis_methods[method](prediction_results, selection)
    raise ValueError(f"Unsupported analysis method: {method}")
