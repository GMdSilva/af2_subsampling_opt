"""
Defines a class to preprocess structural observable data
"""

from typing import Dict, List, Tuple, Union

import numpy as np

from src.utilities.utilities import find_ranges


class DataPreprocessor:
    """
    Class to preprocess structural observable data for downstream purposes.
    """
    def __init__(self, data: Union[List[float], np.array]):
        """
        Args:
        - data (List[float]): The input data list.
        """
        self.data = np.array(data)

    def remove_out(self) -> Dict[str, Union[List[float], List[Tuple[int, int]], float]]:
        """
        Process the data to remove outliers using the IQR method.

        Returns:
        - Dict[str, Union[List[float], List[Tuple[int, int]], float]]:
            A dictionary with the processed data
            (with NaN values for outliers) and the outlier ranges.
        """
        # Calculate Interquartile Range
        quant_1 = np.percentile(self.data, 25)
        quant_3 = np.percentile(self.data, 75)
        iqr = quant_3 - quant_1

        # Set boundaries for outliers
        lower_bound = quant_1 - 1.5 * iqr
        upper_bound = quant_3 + 1.5 * iqr

        # Get boolean masks for valid and outlier data
        is_outlier = (self.data < lower_bound) | (self.data > upper_bound)
        is_valid = ~is_outlier

        # Convert valid data to original values and outliers to NaN
        no_outliers = np.where(is_valid, self.data, np.nan)

        # Find indices of outliers and get their ranges
        outlier_indices = np.where(is_outlier)[0].tolist()
        outlier_ranges = find_ranges(outlier_indices)

        return {
            "no_outliers": no_outliers.tolist(),
            "outlier_ranges": outlier_ranges,
            "lower_bound": lower_bound,
            "upper_bound": upper_bound
        }

    def get_mvavg(self, window_size: int = 10) -> np.array:
        """
        Compute the moving average of a given data set using a specified window size.

        Args:
        - window_size (int, optional): The size of the window for computing the moving average.
            Defaults to 10.

        Returns:
        - np.array: The moving average of the input data.
        """
        kernel = np.ones(window_size) / window_size
        full_convolution = np.convolve(self.data, kernel, mode='full')
        moving_avg = full_convolution[window_size - 1:]
        return moving_avg

    def process_data(self, win_size: int = 10) -> Dict[str, Union[np.array, List[float], float]]:
        """
        Process the provided data by removing outliers and computing its moving average.

        Args:
        - win_size (int, optional): Size of the window for computing the moving avg.

        Returns:
        - dict: A dictionary containing the processed data without outliers and its moving average.
        """
        processed_data = self.remove_out()
        processed_data['moving_avg_no_out'] = self.get_mvavg(win_size)
        return processed_data
