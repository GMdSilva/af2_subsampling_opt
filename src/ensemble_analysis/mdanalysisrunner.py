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

from user_settings.config import FIRST_RESIDUE, REINDEX, PREDICTION_ROOT
from src.utilities.utilities import save_to_pickle


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
        Calculates the Root Mean Square Deviation (RMSD) for given af2_predictions.
            Uses the #1 ranked prediction (by pLDDT) as reference.

        Returns:
        - Dictionary with RMSD analysis results.
        """
        pred_rmsd = rms.RMSD(traj,
                             select='protein and name CA',
                             groupselections=[self.selection])
        pred_rmsd.run()
        return {'residues': pred_rmsd.results.rmsd[:, 1], 'results': pred_rmsd.results.rmsd[:, 3]}

    def calc_rmsd_ref(self, traj) -> Dict[str, np.ndarray]:
        """
        Calculates the RMSD using a reference for given af2_predictions.

        Parameters:
        - reference_path: Path to reference structure.

        Returns:
        - Dictionary with RMSD analysis results using a reference.
        """
        ref = mda.Universe(self.reference)
        ref = ref.select_atoms('protein and name CA')
        pred_rmsd = rms.RMSD(traj,
                             ref,
                             select='protein and name CA',
                             groupselections=[self.selection])
        pred_rmsd.run()
        return {'residues': pred_rmsd.results.rmsd[:, 1], 'results': pred_rmsd.results.rmsd[:, 3]}

    def calc_rmsf(self, traj) -> Dict[str, np.ndarray]:
        """
        Calculates the Root Mean Square Fluctuation (RMSF) for given af2_predictions.

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
        - method: Observable to measure in af2_predictions (e.g. RMSD, distances, etc.)

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

    @staticmethod
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

    @staticmethod
    def load_trajectory(path: str) -> mda.Universe:
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
                        method: str = 'rmsd',
                        label: str = None) -> dict:
        """
        Performs bulk analysis on multiple directories.

        Parameters:
        - method: Observable method for af2_predictions

        Returns:
        - List of dictionaries containing analysis results.
        """
        if self.selection != 'protein and name CA':
            after_resid = self.selection.split("resid", 1)[1]
            # Split the result based on "and name CA" and get the former part
            resid = after_resid.split("and name CA", 1)[0].strip()
            resid = resid.replace(":", "_")
        else:
            resid = ''
        if bulk:
            folders = [d for d in glob(f'{self.path}/{self.prefix}*') if os.path.isdir(d)]
        else:
            folders = [os.path.join(PREDICTION_ROOT,
                                    'results',
                                    'af2_predictions',
                                    f"{self.prefix}_"
                                    f"{trial.split(':')[0]}_"
                                    f"{trial.split(':')[1]}")]
        all_results = []
        for folder in folders:
            print(f"Analyzing {method} of {self.prefix} prediction, "
                  f"parameters {folder.split('_')[-2]}:{folder.split('_')[-1]}, "
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

            if label is not None:
                filename = os.path.join(PREDICTION_ROOT,
                                        'results',
                                        'misc_data',
                                        f'{self.prefix}_'
                                        f'{method}_'
                                        f"{folder.split('_')[-2]}_{folder.split('_')[-1]}_"
                                        f"{resid}_"
                                        f"{label}_"
                                        f"results.pkl")
            else:
                filename = os.path.join(PREDICTION_ROOT,
                                        'results',
                                        'misc_data',
                                        f'{self.prefix}_'
                                        f'{method}_'
                                        f"{folder.split('_')[-2]}_{folder.split('_')[-1]}_"
                                        f"{resid}_"
                                        f"results.pkl")
            save_to_pickle(filename, all_results)

        if bulk:
            filename = os.path.join(PREDICTION_ROOT,
                                    'results',
                                    'misc_data',
                                    f'{self.prefix}_{method}_bulk_results.pkl')
            save_to_pickle(filename, all_results)

        return all_results
