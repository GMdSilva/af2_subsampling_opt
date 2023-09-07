"""
Defines a class to run MDAnalysis
    and measure structural observables in AF2 prediction ensembles.
"""
import glob
import os
from glob import glob
from typing import Dict

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import align, rms

from user_settings.config import FIRST_RESIDUE, REINDEX
from utilities.utilities import save_to_pickle


class MDAnalysisRunner:
    """
    Class to run MDAnalysis and measure structural observables
        in AF2 prediction ensembles.
    """

    def __init__(self, path: str, prefix: str, selection: str = 'protein and name CA', ref_name=""):
        """
        Initialize MDAnalysisRunner instance.

        Parameters:
        - path: Path to folder containing .pdb structures.
        - prefix: Usually the name of the target protein.
        - selection: Atom selection for downstream analyzes.
        - ref_name: Name of the reference.
        """
        self.path = path
        self.prefix = prefix
        self.selection = selection
        self.reference = ref_name

    def calc_rmsd(self, traj) -> Dict[str, np.ndarray]:
        """
        Calculates the Root Mean Square Deviation (RMSD) for given predictions.
            Uses the #1 ranked prediction (by pLDDT) as reference.

        Returns:
        - Dictionary with RMSD analysis results.
        """
        pred_rmsd = rms.RMSD(traj,
                             select='protein and name CA',
                             groupselections=[self.selection])
        pred_rmsd.run()
        return {'residues': pred_rmsd.results.rmsd[:, 1], 'results': pred_rmsd.results.rmsd[:, 3]}

    def calc_rmsd_ref(self, traj, reference_path: str = "") -> Dict[str, np.ndarray]:
        """
        Calculates the RMSD using a reference for given predictions.

        Parameters:
        - reference_path: Path to reference structure.

        Returns:
        - Dictionary with RMSD analysis results using a reference.
        """
        ref = mda.Universe(reference_path)
        ref = ref.select_atoms(self.selection)
        pred_rmsd = rms.RMSD(traj,
                             ref,
                             select='protein and name CA',
                             groupselections=[self.selection])
        pred_rmsd.run()
        return {'residues': pred_rmsd.results.rmsd[:, 1], 'results': pred_rmsd.results.rmsd[:, 3]}

    def calc_rmsf(self, traj) -> Dict[str, np.ndarray]:
        """
        Calculates the Root Mean Square Fluctuation (RMSF) for given predictions.

        Returns:
        - Dictionary with RMSF analysis results.
        """
        average = align.AverageStructure(traj,
                                         traj,
                                         select=self.selection,
                                         ref_frame=0).run()
        ref = average.results.universe
        align.AlignTraj(traj, ref, select=self.selection, in_memory=True).run()
        c_alphas = traj.select_atoms(self.selection)
        pred_rmsf = rms.RMSF(c_alphas).run()
        return {'residues': c_alphas.resids, 'results': pred_rmsf.results.rmsf}

    def analyze_predictions(self, traj: mda.Universe,
                            method: str = 'rmsd') -> Dict[str, np.ndarray]:
        """
        Analyzes a trajectory based on the given configuration.

        Parameters:
        - traj: Trajectory built from observables.
        - method: Observable to measure in predictions (e.g. RMSD, distances, etc.)

        Returns:
        - Dictionary containing analysis results.
        """
        if method == "rmsd":
            return self.calc_rmsd(traj)
        if method == "rmsf":
            return self.calc_rmsf(traj)
        if method == "rmsd_ref":
            return self.calc_rmsd_ref(traj)
        raise ValueError(f"Unsupported analysis method: {method}")

    def extract_trial_name(self, folder: str) -> str:
        """
        Extracts trial name from folder name.

        Parameters:
        - folder: Folder name from which trial name needs to be extracted

        Returns:
        - Extracted trial name
        """
        trial = folder.split('_')
        return f'{trial[-2]}:{trial[-1]}'

    def load_trajectory(self, path: str) -> mda.Universe:
        """
        Load a trajectory from a collection of PDB files using MDAnalysis.

        Parameters:
        - PATH_PATTERN (str): The glob pattern specifying the path to the PDB files.

        Returns:
        - MDAnalysis.Universe: The universe object containing the combined trajectory.
        """

        # Get list of sorted files
        pdb_files = sorted(glob(path + "/*unrelaxed_rank*.pdb"))

        # Use the ChainReader to treat multiple files as a single trajectory
        traj = mda.Universe(pdb_files[0], pdb_files, dt=1)
        if REINDEX:
            for i, residue in enumerate(traj.residues, start=FIRST_RESIDUE):
                residue.resid = i

        return traj

    def process_results(self, bulk: bool = True,
                        trial: str = '256:512',
                        method: str = 'rmsd') -> dict:
        """
        Processes analysis results either in bulk or for a specific trial.

        Parameters:
        - bulk (bool): If True, processes all directories in bulk.
            Otherwise, processes a specific trial.
        - trial (str): Specifies which trial to process if not in bulk mode.
        - method (str): Observable to measure in predictions (e.g. RMSD, distances, etc.)

        Returns:
        - List of dictionaries containing analysis results.
        """
        if bulk:
            all_results = self.bulk_process(method)
        else:
            all_results = []
            print(f"Analyzing {method} of {self.prefix} prediction, "
                  f"parameters: {trial}, "
                  f"with selection: {self.selection}")

            traj = self.load_trajectory(self.path)
            analysis_results = self.analyze_predictions(traj, method=method)

            result_dict = {
                "results": analysis_results['results'],
                'residues': analysis_results['residues'],
                'trial': trial,
                "method": method,
                "selection": self.selection,
                "reference": self.reference,
            }
            all_results.append(result_dict)

            filename = os.path.join('results',
                                    'mdanalysis_results',
                                    f"{self.prefix}_{method}_"
                                    f"{trial.split(':')[0]}_"
                                    f"{trial.split(':')[1]}_"
                                    f"results.pkl")
            save_to_pickle(filename, all_results)

        return all_results

    def bulk_process(self, method: str = 'rmsd') -> list:
        """
        Performs bulk analysis on multiple directories.

        Parameters:
        - method: Observable method for predictions

        Returns:
        - List of dictionaries containing analysis results.
        """
        folders = [d for d in glob(f'{self.path}/{self.prefix}*') if os.path.isdir(d)]
        all_results = []

        for folder in folders:
            print(f"Analyzing {method} of {self.prefix} prediction, "
                  f"max_seq {folder.split('_')[-2]}, "
                  f"with selection: {self.selection}")
            trial = self.extract_trial_name(folder)
            traj = self.load_trajectory(folder)
            analysis_results = self.analyze_predictions(traj, method=method)

            result_dict = {
                "results": analysis_results['results'],
                'residues': analysis_results['residues'],
                'trial': trial,
                "method": method,
                "selection": self.selection,
                "reference": self.reference,
            }
            all_results.append(result_dict)

        filename = os.path.join('results',
                                'mdanalysis_results',
                                f'{self.prefix}_{method}_results.pkl')
        save_to_pickle(filename, all_results)

        return all_results
