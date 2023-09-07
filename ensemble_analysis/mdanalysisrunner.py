"""
Defines a class to run MDAnalysis
    and measure structural observables in AF2 prediction ensembles.
"""

from typing import Dict

import MDAnalysis
from MDAnalysis.analysis import rms, align
import numpy as np

from utilities.utilities import load_trajectory


class MDAnalysisRunner:
    """
    Class to run MDAnalysis and measure structural observables
        in AF2 prediction ensembles.
    """
    def __init__(self, path: str, selection: str = 'protein and name CA'):
        """
        Parameters:
        - path: Path to folder containing .pdb structures.
        - selection: Atom selection for downstream analyzes.
        """
        self.path = path
        self.selection = selection
        self.prediction_results = load_trajectory(path)

    def calc_rmsd(self) -> Dict[str, np.ndarray]:
        """
        Calculates the Root Mean Square Deviation (RMSD) for given predictions.
            Uses the #1 ranked prediction (by pLDDT) as reference.

        Returns:
        - Dictionary with RMSD analysis results.
        """
        pred_rmsd = rms.RMSD(self.prediction_results,
                             select='protein and name CA',
                             groupselections=[self.selection])
        pred_rmsd.run()
        return {'residues': pred_rmsd.results.rmsd[:, 1], 'results': pred_rmsd.results.rmsd[:, 3]}

    def calc_rmsd_ref(self, reference_path: str) -> Dict[str, np.ndarray]:
        """
        Calculates the RMSD using a reference for given predictions.

        Parameters:
        - reference_path: Path to reference structure.

        Returns:
        - Dictionary with RMSD analysis results using a reference.
        """
        ref = MDAnalysis.Universe(reference_path)
        ref = ref.select_atoms(self.selection)
        pred_rmsd = rms.RMSD(self.prediction_results,
                             ref,
                             select='protein and name CA',
                             groupselections=[self.selection])
        pred_rmsd.run()
        return {'residues': pred_rmsd.results.rmsd[:, 1], 'results': pred_rmsd.results.rmsd[:, 3]}

    def calc_rmsf(self) -> Dict[str, np.ndarray]:
        """
        Calculates the Root Mean Square Fluctuation (RMSF) for given predictions.

        Returns:
        - Dictionary with RMSF analysis results.
        """
        average = align.AverageStructure(self.prediction_results,
                                         self.prediction_results,
                                         select=self.selection,
                                         ref_frame=0).run()
        ref = average.results.universe
        align.AlignTraj(self.prediction_results, ref, select=self.selection, in_memory=True).run()
        c_alphas = self.prediction_results.select_atoms(self.selection)
        pred_rmsf = rms.RMSF(c_alphas).run()
        return {'residues': c_alphas.resids, 'results': pred_rmsf.results.rmsf}

    def analyze_predictions(self, method: str = 'rmsd') -> Dict[str, np.ndarray]:
        """
        Analyzes a trajectory based on the given configuration.

        Parameters:
        - method: observable to measure in predictions (e.g. RMSD, distances, etc.)

        Returns:
        - Dictionary containing analysis results.
        """
        analysis_methods = {
            'rmsd': self.calc_rmsd,
            'rmsd_ref': self.calc_rmsd_ref,
            'rmsf': self.calc_rmsf
        }

        if method in analysis_methods:
            return analysis_methods[method]()
        raise ValueError(f"Unsupported analysis method: {method}")
