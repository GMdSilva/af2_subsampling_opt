""" Misc. scripts for different purposes """
from typing import List, Tuple, Any
from pathlib import Path
import pickle
import glob
import MDAnalysis as mda

from user_settings.config import REINDEX, FIRST_RESIDUE


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


def load_trajectory(path: str) -> mda.Universe:
    """
    Load a trajectory from a collection of PDB files using MDAnalysis.

    Parameters:
    - PATH_PATTERN (str): The glob pattern specifying the path to the PDB files.

    Returns:
    - MDAnalysis.Universe: The universe object containing the combined trajectory.
    """

    # Get list of sorted files
    pdb_files = sorted(glob.glob(path + "/*unrelaxed_rank*.pdb"))

    # Use the ChainReader to treat multiple files as a single trajectory
    traj = mda.Universe(pdb_files[0], pdb_files, dt=1)
    if REINDEX:
        for i, residue in enumerate(traj.residues, start=FIRST_RESIDUE):
            residue.resid = i

    return traj


def tuple_range_to_string(range_tuple: Tuple[int, int]) -> str:
    """
    Convert a tuple representing a range to a string in the format "start-end".

    Args:
        range_tuple (Tuple[int, int]): A tuple representing a range as (start, end).

    Returns:
        str: String representation of the range in the format "start-end".
    """
    start, end = range_tuple
    return f"{start}-{end}"


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
