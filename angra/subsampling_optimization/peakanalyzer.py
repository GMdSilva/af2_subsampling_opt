"""
Defines a class that analyzes modes/peaks
    obtained from each subsampling parameter set.
"""

import os
from typing import Dict, List, Tuple, Union

import numpy as np
from scipy.signal import find_peaks
from sklearn.neighbors import KernelDensity

from angra.utilities.utilities import save_to_pickle
from user_settings.config import PREDICTION_ROOT


class PeakAnalyzer:
    """
    Class that analyzes modes/peaks
        obtained from each subsampling parameter set.
    """

    def __init__(self, prefix: str, trial: str):
        """
        Args:
            - prefix (str): Usually the name of the target protein.
            - trial (str): Specific parameter set being tested.
        """
        self.prefix = prefix
        self.trial = trial

    @staticmethod
    def get_silverman_bandwidth(data: Union[List[float], np.ndarray]) -> float:
        """
        Calculate the bandwidth using Silverman's Rule of Thumb.

        Args:
            - data (Union[List[float], np.ndarray]): Dataset for which to calculate the bandwidth.

        Returns:
            float: The calculated bandwidth.
        """
        sample_size = len(data)
        std_dev = np.std(data, ddof=1)
        silverman_factor = (4 * (std_dev ** 5) / (3 * sample_size))
        return silverman_factor ** 0.2

    @staticmethod
    def estimate_density(data: Union[List[float], np.ndarray],
                         bandwidth: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Estimate the density using KDE.

        Args:
            - data (Union[List[float], np.ndarray]): Input data for density estimation.
            - bandwidth (float): Bandwidth for Kernel Density Estimation.

        Returns:
            - Tuple[np.ndarray, np.ndarray]: x values and their corresponding density estimates.
        """
        kde = KernelDensity(bandwidth=bandwidth, kernel='gaussian')
        kde.fit(np.array(data).reshape(-1, 1))
        x_values = np.linspace(min(data), max(data), 1000)
        density = np.exp(kde.score_samples(x_values.reshape(-1, 1)))
        return x_values, density

    @staticmethod
    def find_modalities(density: Union[List[float], np.ndarray],
                        x_values: Union[List[float], np.ndarray],
                        qtl: float = 30) -> Tuple[int, np.ndarray, np.ndarray, np.ndarray, float]:
        """
        Find the modalities using a quantile-based threshold.

        Args:
            - density (Union[List[float], np.ndarray]): Density values corresponding to x values.
            - x_values (Union[List[float], np.ndarray]): x values where the density is estimated.
            - qtl (float, optional): Quantile used for setting the threshold. Defaults to 30.

        Returns:
            Tuple[int, np.ndarray, np.ndarray, np.ndarray, float]:
            - num_modes (int): Number of modalities.
            - populations (np.ndarray): Heights of the peaks.
            - modes (np.ndarray): x values corresponding to the peaks.
            - peaks (np.ndarray): Indices of the peaks.
            - threshold (float): Threshold used for peak detection.
        """
        threshold = np.percentile(density, qtl)
        peaks, properties = find_peaks(density, height=threshold)

        return len(peaks), properties['peak_heights'], x_values[peaks], peaks, threshold

    @staticmethod
    def sort_peaks_by_population(mode_data: Dict[str, Union[np.ndarray, List]]) \
            -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Sort peaks based on their populations.

        Args:
            - mode_data (Dict[str, Union[np.ndarray, List]]): Dictionary containing:
                - 'populations': Array of peak populations.
                - 'peaks': Indices of peaks in the distribution.
                - 'distribution': Array of x-values where the density is estimated.
                - 'modes': x values corresponding to the peaks.

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray]:
            - sorted_pops (np.ndarray): Sorted populations in ascending order.
            - sorted_peaks (np.ndarray): x values of sorted peaks.
            - sorted_modes (np.ndarray): x values of sorted modes.
        """
        pops = mode_data['populations']
        peaks = mode_data['peaks']
        x_values = mode_data['distribution']
        modes = mode_data['modes']
        sorted_indices = np.argsort(pops)
        return pops[sorted_indices], x_values[peaks][sorted_indices], modes[sorted_indices]

    def analyze_modes(self,
                      data: Union[List[float], np.ndarray],
                      bandwidth: float = 0.1,
                      quantile: float = 30) ->\
            Dict[str, Union[int, float, List[float]]]:
        """
        Analyze the modality of a dataset using KDE and a quantile-based threshold.

        Args:
            - data (Union[List[float], np.ndarray]): Dataset for which to estimate modality.
            - rmsd_range (Union[List[float], np.ndarray]): Range for RMSD.
            - bandwidth (float, optional): Bandwidth for KDE. Defaults to 0.1.
            - quantile (float, optional): Percentile for peak detection. Defaults to 30.

        Returns:
            dict: With keys 'num_modes', 'populations',
                'modes', 'peaks', 'distribution', and 'trial'.
        """

        x_values, density = self.estimate_density(data, bandwidth)
        num_modes, pops, modes, peaks, threshold = self.find_modalities(density, x_values, quantile)
        return {
            "density": density,
            "x_values": x_values,
            "threshold": threshold,
            "num_modes": num_modes,
            "populations": pops,
            "modes": modes,
            'peaks': peaks,
            'distribution': x_values,
            'trial': self.trial
        }

    def two_state_analysis(self,
                           mode_data: Dict[str, Union[int, float, List[float]]],
                           residue_range: str) -> Dict[str, Union[int, float, List[float]]]:
        """
        Analyze mode data for a system with two states.

        Args:
            - mode_data (dict): Dictionary containing mode information.
            - rmsd_range (str): Residue range corresponding to the analyzed data

        Returns:
            - dict: Analyzed results.
        """
        sorted_pops, sorted_peaks, sorted_modes = self.sort_peaks_by_population(mode_data)
        ground_state_population = sorted_pops[-1]
        alt1_state_population = sorted_pops[-2]
        ground_state_peak = sorted_peaks[-1]
        alt1_state_peak = sorted_peaks[-2]
        distance_between_peaks = alt1_state_peak - ground_state_peak
        return {
            'num_modes': mode_data['num_modes'],
            'populations': mode_data['populations'],
            'modes': sorted_modes,
            'trial': mode_data['trial'],
            'ground_pop': ground_state_population,
            'alt1_pop': alt1_state_population,
            'ground_peak': ground_state_peak,
            'alt1_peak': alt1_state_peak,
            'dist_ground_alt1': distance_between_peaks,
            'distribution': mode_data['distribution'],
            'residue_range': residue_range
        }

    def mode_sanity_checks(self, modality_results: dict, rmsd_range: str) -> bool:
        """
        Determines if the modality_results pass certain conditions.

        Conditions:
        1. The number of modes should be between 2 and 3 (inclusive).
        2. For all modes, the population should not exceed the population of the first mode.
        3. If any mode's distance exceeds the threshold, it is considered "passing".

        Args:
            - modality_results (dict): Dictionary containing 'num_modes',
                'populations', 'distances', and 'thresh'.
            - rmsd_range (str): range of the residues used to calc. RMSD.

        Returns:
            - bool: True if the results pass the conditions, otherwise False.
        """
        num_modes = modality_results['num_modes']
        ground_peak = modality_results['ground_peak']
        alt1_peak = modality_results['alt1_peak']
        max_value = modality_results['distribution'].max()
        if max_value > 17:
            print(f"Peak range {rmsd_range} with Alt1 RMSD of {alt1_peak:.3g} A "
                  f"for parameter {self.trial} rejected: "
                  f"RMSD distribution too wide, likely unfolded. Max RMSD value: ({max_value}).")
            return False
        if not 1 <= num_modes <= 3:
            print(f"Peak range {rmsd_range} with Alt1 RMSD of {alt1_peak:.3g} A "
                  f"for parameter {self.trial} rejected: "
                  f"too many alternative conformations modes ({num_modes}).")
            return False
        if alt1_peak > 12:
            print(f"Peak range {rmsd_range} with Alt1 RMSD of {alt1_peak:.3g} A "
                  f"for parameter {self.trial} rejected: RMSD > 12, likely unfolded.")
            return False
        if ground_peak > alt1_peak:
            print(f"Peak range {rmsd_range} with Alt1 RMSD of {alt1_peak:.3g} A "
                  f"for parameter {self.trial} rejected: "
                  f"ground conformation has higher RMSD than alt1 conformation.")
            return False
        # if alt1_peak - ground_peak < 1.5:
        #     print(f"Peak range {rmsd_range} with Alt1 RMSD of {alt1_peak:.3g} A "
        #           f"for parameter {self.trial} rejected: "
        #           f"Ground conformation too close to Alt1 conformation.")
        #     return False

        print(f"Peak range {rmsd_range} with Alt1 RMSD of {alt1_peak:.3g} A "
              f"for parameter {self.trial} accepted ")

        filename = f"{self.prefix}_accepted_" \
                   f"{self.trial.split(':')[0]}_" \
                   f"{self.trial.split(':')[1]}_" \
                   f"range_{rmsd_range[0]}_" \
                   f"{rmsd_range[1]}.pkl"

        peaks_path = os.path.join(PREDICTION_ROOT,
                                  "results",
                                  "optimization_results",
                                  "peaks",
                                  filename)

        save_to_pickle(peaks_path, modality_results)
        return True
