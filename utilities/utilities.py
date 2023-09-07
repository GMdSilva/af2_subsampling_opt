""" Misc. scripts for different purposes """
import os
import pickle
from pathlib import Path
from typing import Any, List, Tuple


def load_from_pickle(filename: Path) -> Any:
    """
    Load data from a pickle file.

    Args:
        filename (Path): The path to the pickle file to load.

    Returns:
        Any: The deserialized object from the pickle file.
    """
    with open(f"{filename}", 'rb') as file:
        return pickle.load(file)


def save_to_pickle(filename: Path, data: Any) -> None:
    """
    Save data to a pickle file.

    Args:
        filename (Path): The path where to save the pickle file.
        data (Any): The data object to serialize and save.
    """
    with open(f"{filename}", 'wb') as file:
        pickle.dump(data, file)


def tuple_range_to_string(range_tuple: Tuple[int, int]) -> str:
    """
    Convert a tuple representing a range to a string in the format "start-end".

    Args:
        range_tuple (Tuple[int, int]): A tuple representing a range as (start, end).

    Returns:
        str: String representation of the range in the format "start-end".
    """
    start, end = range_tuple
    return f"{start}:{end}"


def ranges_overlap(range_a: Tuple[int, int],
                   range_b: Tuple[int, int],
                   wiggle: int = 0) -> bool:
    """
    Check if two given ranges overlap considering an additional wiggle room.

    Args:
        range_a (Tuple[int, int]): A tuple representing the first range as (start, end).
        range_b (Tuple[int, int]): A tuple representing the second range as (start, end).
        wiggle (int, optional): Wiggle room to consider when checking overlap. Defaults to 0.

    Returns:
        bool: True if ranges overlap with the wiggle room, otherwise False.
    """
    start_a, end_a = range_a
    start_b, end_b = range_b
    return end_a + wiggle >= start_b and end_b + wiggle >= start_a


def find_common_ranges_with_wiggle(all_ranges: List[List[Tuple[int, int]]],
                                   wiggle: int = 5) -> List[Tuple[int, int]]:
    """
    Find common ranges between multiple sets of ranges considering a wiggle room.

    Args:
        all_ranges (List[List[Tuple[int, int]]]): List of list of ranges, each as (start, end).
        wiggle (int, optional): The wiggle room to consider when checking. Defaults to 5.

    Returns:
        List[Tuple[int, int]]: List of common ranges.
    """

    common_ranges = set(all_ranges[0])

    for curr_ranges in all_ranges[1:]:
        overlaps = set()
        for curr_range in curr_ranges:
            overlaps.update({r for r in common_ranges if ranges_overlap(r, curr_range, wiggle)})
        common_ranges &= overlaps

    return list(common_ranges)


def find_ranges(indices_list: List[int]) -> List[Tuple[int, int]]:
    """
    Convert a list of indices into consecutive ranges.

    Args:
    - indices_list (List[int]): List of indices.

    Returns:
    - List[Tuple[int, int]]: A list of tuples representing the start and end of consecutive ranges.
    """
    ranges = []

    if not indices_list:
        return ranges

    start_idx = indices_list[0]
    end_idx = indices_list[0]

    for idx in indices_list[1:]:
        if idx == end_idx + 1:
            end_idx = idx
        else:
            ranges.append((start_idx, end_idx))
            start_idx = end_idx = idx

    ranges.append((start_idx, end_idx))
    return ranges


def rename_files_recursively(directory_path: str, old_prefix: str, new_prefix: str) -> None:
    """
    Rename files in a directory and its subdirectories from a specified prefix to a new prefix,
    specifically for files with names following the pattern: <prefix>_unrelaxed_rank_*.pdb.

    Parameters:
    - directory_path (str): Path to the root directory containing the .pdb files.
    - old_prefix (str): The old prefix in the filenames to be replaced.
    - new_prefix (str): The new prefix to replace the old one.

    Returns:
    - None: The function prints out the renamed files.
    """

    # Use os.walk to iterate over each directory and its subdirectories
    for dirpath, dirnames, filenames in os.walk(directory_path):
        for filename in filenames:
            # Check if the file has a .pdb extension
            if filename.endswith(".pdb"):
                # Construct the old pattern using the old prefix
                old_pattern = f"{old_prefix}_unrelaxed_rank_"

                # Check if the filename matches the old pattern
                if old_pattern in filename:
                    # Construct the new pattern using the new prefix
                    new_pattern = f"{new_prefix}_unrelaxed_rank_"

                    # Replace the old pattern with the new pattern in the filename
                    new_filename = filename.replace(old_pattern, new_pattern)

                    # Construct the full path for the old and new filenames
                    old_file_path = os.path.join(dirpath, filename)
                    new_file_path = os.path.join(dirpath, new_filename)

                    # Rename the file
                    os.rename(old_file_path, new_file_path)

                    # Print the renaming action for confirmation
                    print(f"Renamed {old_file_path} to {new_file_path}")
