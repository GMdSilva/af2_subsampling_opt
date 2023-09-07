"""
Defines a class to find top n peaks and calculate scores for them.
"""

from typing import List, Tuple, Dict

import numpy as np
from scipy.signal import find_peaks

from ensemble_analysis.datapreprocessor import DataPreprocessor


class PeakDetector:
    """
    Class to find top n residue ranges and calculate scores for them.
    """
    def __init__(self, result: Dict[str, List[float]]):
        """
        Args:
        - result (Dict[str, List[float]]):
            Dictionary containing 'results' key that has data for processing.
        """
        self.result = result

    def get_peaks(self, peak_thresh: float) -> List[Dict[str, float]]:
        """
        Detect and sort peaks in the smoothed analysis results.

        Parameters:
        - peak_thresh (float): Threshold for peak detection.

        Returns:
        - List[Dict[str, float]]: List of dictionaries with peak properties.
        """
        peaks, properties = find_peaks(self.result['results_smooth'], width=(peak_thresh, None))
        peak_heights = [self.result['results_smooth'][idx] for idx in peaks]
        sorted_indices = np.argsort(peak_heights)[::-1]
        peak_results = [
            {
                'height': self.result['results_smooth'][peaks[idx]],
                'left_base': properties['left_bases'][idx],
                'right_base': properties['right_bases'][idx]
            }
            for idx in sorted_indices
        ]
        return peak_results

    @staticmethod
    def trapezoid_area_from_points(x_values: List[float], y_values: List[float]) -> float:
        """
        Estimate the area under the curve using the trapezoid rule.

        Parameters:
        - x_values (List[float]): The x values.
        - y_values (List[float]): The y values corresponding to each x value.

        Returns:
        - float: Estimated area under the curve.
        """
        return sum(0.5 * (x1 - x0) * (y1 + y0)
                   for x0, x1, y0, y1 in zip(x_values[:-1],
                                             x_values[1:],
                                             y_values[:-1],
                                             y_values[1:]))

    def normalized_area_in_interval(self,
                                    x_values: List[float],
                                    y_values: List[float],
                                    int_start: float,
                                    int_end: float) -> float:
        """
        Calculate the normalized area under the curve
            in the interval [a, b] using the trapezoid rule.

        Parameters:
        - x_values (List[float]): The x values.
        - y_values (List[float]): The y values corresponding to each x.
        - int_start (float): The start of the interval.
        - int_end (float): The end of the interval.

        Returns:
        - float: Normalized area in the interval [int_start, int_end].
        """
        filtered_x = [xi for xi in x_values if int_start <= xi <= int_end]
        filtered_y = [y_values[i] for i, xi in enumerate(x_values) if int_start <= xi <= int_end]
        raw_area = self.trapezoid_area_from_points(filtered_x, filtered_y)
        interval_length = int_end - int_start
        normalized_area = raw_area / interval_length
        return normalized_area

    def sort_ranges_by_auc(self, ranges: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
        """
        Sort ranges based on the normalized area under the curve in descending order.

        Parameters:
        - ranges (List[Tuple[int, int]]): List of tuple ranges.
        - result (Dict[str, List[float]]): Dictionary containing 'residues' and 'results' lists.

        Returns:
        - List[Tuple[int, int]]: Ranges sorted by the normalized area in descending order.
        """
        interval_areas = [
            self.normalized_area_in_interval(self.result['residues'],
                                             self.result['results'], start, end)
            for start, end in ranges
        ]
        sorted_indices = np.argsort(-np.array(interval_areas))
        return [ranges[i] for i in sorted_indices]

    @staticmethod
    def refactor_peaks(peak_results: List[Dict[str, int]]) -> List[Tuple[int, int]]:
        """
        Gets peak_ranges from dict containing left and right bases

        Args:
        - peak_results (List[Dict[str, int]]): List of dictionaries with peak results.

        Returns:
        - List[Tuple[int, int]]: Combined list of peak ranges.
        """
        peak_ranges = [(int(peak['left_base']),
                        int(peak['right_base']))
                       for peak in peak_results]
        return peak_ranges

    def detect_peaks(self,
                     window_size: int = 5,
                     peak_thresh: int = 3,
                     peak_n: int = 5,
                     smooth: bool = False) -> List[Tuple[int, int]]:
        """
        Detects peaks from given observable analysis results

        Args:
        - window_size (int): Size of the window for moving average smoothing.
        - peak_thresh (int): Threshold for peak detection.
        - peak_n (int): Minimum number of residues required to form a peak.
        - smooth (bool): If true, uses moving average

        Returns:
        - list: List of largest overlapping ranges representing peaks.
        """
        data_processor = DataPreprocessor(data=self.result['results'])
        if smooth:
            data = data_processor.get_mvavg(window_size=window_size)
        else:
            data = self.result['results']
        self.result['results_smooth'] = data
        peaks = self.get_peaks(peak_thresh)
        peak_ranges = self.refactor_peaks(peaks)
        largest_ranges = peak_ranges
        top_n_peak_ranges = self.sort_ranges_by_auc(largest_ranges)[:peak_n]
        return top_n_peak_ranges
