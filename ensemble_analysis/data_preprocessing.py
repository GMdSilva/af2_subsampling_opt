"""Preprocesses structural observable data for downstream purposes"""

from typing import List, Dict, Union, Tuple

import numpy as np
import matplotlib.pyplot as plt

from utilities.utilities import find_ranges


def remove_out(data: List[float]) -> Dict[str, Union[List[float], List[Tuple[int, int]], float]]:
    """
    Process the data to remove outliers using the IQR method.

    Args:
    - data (List[float]): The input data list.

    Returns:
    - Dict[str, Union[List[float], List[Tuple[int, int]], float]]:
        A dictionary with the processed data
        (with NaN values for outliers) and the outlier ranges.
    """

    # Convert data to numpy array for efficient operations
    data_array = np.array(data)

    # Calculate Interquartile Range
    quant_1 = np.percentile(data_array, 25)
    quant_3 = np.percentile(data_array, 75)
    iqr = quant_3 - quant_1

    # Set boundaries for outliers
    lower_bound = quant_1 - 1.5 * iqr
    upper_bound = quant_3 + 1.5 * iqr

    # Get boolean masks for valid and outlier data
    is_outlier = (data_array < lower_bound) | (data_array > upper_bound)
    is_valid = ~is_outlier

    # Convert valid data to original values and outliers to NaN
    no_outliers = np.where(is_valid, data_array, np.nan)

    # Find indices of outliers and get their ranges
    outlier_indices = np.where(is_outlier)[0].tolist()
    outlier_ranges = find_ranges(outlier_indices)

    return {
        "no_outliers": no_outliers.tolist(),
        "outlier_ranges": outlier_ranges,
        "lower_bound": lower_bound,
        "upper_bound": upper_bound
    }


def get_mvavg(data: Union[List[float], np.array], window_size: int = 10) -> np.array:
    """
    Compute the moving average of a given data set using a specified window size.

    Args:
    - data (Union[List[float], np.array]): The input data for which the moving average
        is to be computed.
    - window_size (int, optional): The size of the window for computing the moving average.
        Defaults to 10.

    Returns:
    - np.array: The moving average of the input data.
    """

    # Create an array of ones of size window_size and normalize it to sum up to 1.
    # This will be our filter for convolution.
    kernel = np.ones(window_size) / window_size

    # Use convolution to compute the moving average.
    # mode='full' gives the full convolution result.
    # To make the result the same size as the input data, we'll later slice the result.
    full_convolution = np.convolve(data, kernel, mode='full')

    # Slice the convolution result to make it the same size as the input data.
    # For a window_size of 10, the first 9 points don't have 10 points before them.
    moving_avg = full_convolution[window_size - 1:]

    return moving_avg


def process_data(data: Union[List[float], np.array], window_size: int = 10) \
        -> Dict[str, Union[np.array, List[float], float]]:
    """
    Process the provided data by removing outliers and computing its moving average.

    Args:
    - data (list or np.array): The input data to be processed.
    - window_size (int, optional): Size of the window for computing the moving avg. Defaults to 10.

    Returns:
    - dict: A dictionary containing the processed data without outliers and its moving average.
    """

    # Remove outliers from the data
    processed_data = remove_out(data)

    # Compute the moving average for the data without outliers
    processed_data['moving_avg_no_out'] = get_mvavg(processed_data['no_outliers'], window_size)
    plt.plot(processed_data['moving_avg_no_out'])
    plt.show()
    return processed_data
