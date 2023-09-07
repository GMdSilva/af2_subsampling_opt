"""
Defines a class for running subsampled AF2 with colabfold_batch.
"""

import os
from user_settings.config import COLABFOLDBATCH_PATH


class AF2Runner:
    """
    Class for running subsampled AF2 with colabfold_batch.
    """
    def __init__(self,
                 prefix: str,
                 msa: str,
                 trials: int = 5,
                 seeds: int = 32,):
        """
        Parameters:
        - prefix: Usually the name of the target protein.
        - msa: Path to text file containing MSA for prediction.
        - trials: Number of trials to run optimization for.
        - seeds: Number of seeds per prediction.
        - kind: Wild-type (wt) vs. Mutant (e.g. v381v)
        """
        self.prefix = prefix
        self.msa = msa
        self.trials = trials
        self.seeds = seeds

    @staticmethod
    def estimate_starting_parameters(msa: str) -> int:
        """
        Estimate starting parameters based on MSA depth.

        Args:
        - msa (str): Input multi-sequence alignment.

        Returns:
        - int: Starting exponent for sequence count.
        """
        msa_depth = int(msa.split("_")[-2])
        if msa_depth <= 2000:
            return 2
        if 2000 < msa_depth <= 20000:
            return 3
        return 4

    def build_parameter_range(self) -> list:
        """
        Construct a range of parameters for AF2 based on trials and MSA.

        Returns:
        - list: List of parameter sets.
        """
        starting_max_seq_exp = self.estimate_starting_parameters(self.msa)
        parameters = [
            f"{2 ** exp}:{2 ** (exp + 1)}"
            for exp in range(starting_max_seq_exp, self.trials + starting_max_seq_exp)
        ]
        return parameters

    def run_af2(self, parameters: list):
        """
        Execute the colabfold_batch process with the given parameters.

        Args:
        - parameters (list): List of parameter sets.
        """
        result_path = os.path.join('../results', 'predictions')
        for p_set in parameters:
            max_seq, extra_seq = p_set.split(':')
            command = (
                f"colabfold_batch --num-seeds {self.seeds} --max-msa {p_set} --use-dropout "
                f"{self.msa} {os.path.join(result_path, f'{self.prefix}_{max_seq}_{extra_seq}')}"
            )
            os.system(command)

    def make_predictions(self):
        """
        Generate predictions.
        """
        if not os.path.isfile(COLABFOLDBATCH_PATH):
            print('colabfoldbatch not found, not making predictions')
            return

        parameters = self.build_parameter_range()
        self.run_af2(parameters)

    def run_subsampled_af2(self):
        """
        Entry point function to generate subsampled AF2 predictions.
        """
        self.make_predictions()
