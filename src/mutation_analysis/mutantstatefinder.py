from typing import Dict, Tuple, List
import glob

from src.subsampling_optimization.statefinder import StateFinder

class MutantStateFinder(StateFinder):
    """
    Class for analyzing mutations and measuring accuracy at predicting the effects
    of mutations in the ground or alternative states.
    """

    def __init__(self, prefix: str):
        """
        Initializes the MutationAnalyzer class.

        Args:
            prefix (str): Name of the protein being studied.
        """
        self.prefix = prefix
        self.selection: str = 'protein and name CA'
        self.trial: str = None
        self.all_prefixes: list = self._get_prefixes()

    def get_refs_and_compare_muts(self,
                             optimization_results: Dict[str, int]) -> Tuple[List[str], List[str]]:
        """
        Finds PDB files based on optimization results.

        Args:
            optimization_results (Dict[str, int]): A dictionary containing indexes
                                                  for ground and alt1 references.

        Returns:
            Tuple[List[str], List[str]]: A tuple containing two lists.
                                         The first list contains matching paths for ground,
                                         and the second list contains matching paths for alt1.
        """
        indexes = self.get_reps(optimization_results)

        ground_ref_path = self._construct_path(indexes['ground_ref_index'])
        alt1_ref_path = self._construct_path(indexes['alt1_ref_index'])

        ground_files = glob.glob(ground_ref_path)

        alt1_files = glob.glob(alt1_ref_path)

        labeled_files = list(zip(ground_files,
                                 ["ground"] * len(ground_files))) + \
                        list(zip(alt1_files, ["alt1"] * len(alt1_files)))

        rmsds = {}
        for file, label in labeled_files:
            self._save_rep_states(file, label)
            rmsd = self.get_rmsd_vs_refs(file, label)
            rmsds[label] = rmsd

        self._plot_results(rmsds)
        return rmsds

