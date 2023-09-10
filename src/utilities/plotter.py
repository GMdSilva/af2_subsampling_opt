"""
Defines a class for plotting
"""

import os
from typing import Dict, List, Any, Union, Tuple
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot as plt

from user_settings.config import PREDICTION_ROOT


class Plotter:
    """
    Class for plotting
    """
    def __init__(self, prefix: str, data: Union[List, Dict]):
        """
        Initialize the SeabornPlotter with a dataset.

        Args:
        - data (Union(List, Dict)): List or dict containing the data to plot.
        - prefix (str): Usually, the name of the target protein
        """
        self.data = data
        self.prefix = prefix
        plt.close()

        plt.rcParams['xtick.labelsize'] = 16  # x ticks
        plt.rcParams['ytick.labelsize'] = 16  # y ticks
        plt.rcParams['axes.labelsize'] = 18  # x and y labels
        plt.rcParams['axes.titlesize'] = 20  # title of each subplot
        plt.rcParams['figure.titlesize'] = 20  # overall figure title
        plt.rcParams['legend.fontsize'] = 16  # figure legend fontsize

    def plot_rmsf_peaks(self, peak_range: list) -> None:
        """
        Plots peaks for the given results.

        Parameters:
        - peak_range: A list of tuple ranges indicating peaks.
        """
        plt.close()
        colors = ['red',
                  'blue',
                  'green',
                  'yellow',
                  'purple',
                  'orange',
                  'teal',
                  'fuchsia',
                  'pink']

        plt.plot(self.data['residues'], self.data['results'],
                 marker=None,
                 label=self.data['trial'],
                 color='grey')
        plt.xlabel('Residue #')
        plt.ylabel('C-Alpha RMSF (A)')
        plt.title(self.prefix)
        plt.legend()

        for i, (start, end) in enumerate(peak_range):
            color = colors[i % len(colors)]
            plt.plot(self.data['residues'][start:end],
                     self.data['results'][start:end],
                     marker='*',
                     color=color,
                     linestyle='-',
                     markersize=plt.rcParams['lines.markersize'] / 2)
        plt.tight_layout()
        save_path = os.path.join(PREDICTION_ROOT,
                                 'results',
                                 'plots',
                                 f"{self.prefix}_rmsf_trial_"
                                 f"{self.data['trial'].split(':')[0]}_"
                                 f"{self.data['trial'].split(':')[1]}.png")
        plt.savefig(save_path)

    def plot_kde_with_modes(self, peak_range: Tuple[int, int]) -> None:
        """
        Plots Kernel Density Estimation (KDE) with mode peaks.

        This function generates a KDE plot based on the provided data,
        highlights the mode peaks in red, and saves the plot to the specified path.

        Parameters:
        - peak_range (tuple): A tuple of two integers representing the range of peaks.

        Returns:
        None. The function saves the generated plot to the specified path.
        """
        plt.plot(self.data["x_values"],
                 self.data["density"],
                 label=self.data["trial"])
        plt.scatter(self.data["x_values"][self.data["peaks"]],
                    self.data["density"][self.data["peaks"]],
                    color='red')
        plt.xlabel('Value (A)')
        plt.ylabel('Density')
        plt.title(f'Residue Range {str(peak_range)}')
        plt.axhline(y=self.data['threshold'],
                    color='green',
                    linestyle='--')
        plt.legend()
        plt.tight_layout()
        save_path = os.path.join(PREDICTION_ROOT,
                                 'results',
                                 'plots',
                                 f"{self.prefix}_modes_trial_"
                                 f"{self.data['trial'].split(':')[0]}_"
                                 f"{self.data['trial'].split(':')[1]}.png"
                                 f"_ranges_{peak_range[0]}"
                                 f"_{peak_range[1]}.png")
        plt.savefig(save_path)

    def plot_mut_analysis(self, x_values: Dict, mut_list: Dict) -> None:
        """
        Generate a bar plot based on mutation analysis.

        Parameters:
        - x_values (dict): A dictionary containing the key-value pairs for the x-axis.
        - mut_list (dict): A dictionary containing mutation details.

        Returns:
        None. A bar plot is saved in the designated path.
        """

        ground_pop_diffs = [d[self.prefix][list(x_values.keys())[0]]
                            for d in self.data for self.prefix in d]
        labels = [entry['label'] for entry in mut_list.values()]
        effects_ground = [entry['effect']['ground_pop']
                          for entry in mut_list.values()]
        effects_alt1 = [entry['effect']['alt1_pop']
                        for entry in mut_list.values()]

        if x_values == 'ground_pop_diff' or 'ground_pop_test':
            colors = ['red' if effect == '-'
                      else ('grey' if effect == "ref"
                            else 'blue') for effect in effects_ground]
        else:
            colors = ['red' if effect == '-'
                      else ('grey' if effect == "ref"
                            else 'blue') for effect in effects_alt1]

        # Plotting
        fig, ax = plt.subplots()
        ax.barh(labels,
                ground_pop_diffs,
                color=colors)
        ax.set_xlabel(list(x_values.values())[0])
        ax.set_ylabel('')
        ax.set_title('')

        plt.tight_layout()
        save_path = os.path.join(PREDICTION_ROOT,
                                 'results',
                                 'plots',
                                 f"{self.prefix}_{list(x_values.keys())[0]}.png")
        plt.savefig(save_path)


    def plot_multiple_scatter(self, data_list, labels=None, annotations=None):
        """
        Creates a separate scatter plot for each dataset, organized in an automatically determined grid.

        Args:
        - data_list (list of tuples): A list where each element is a tuple containing two lists (x_values, y_values).
        - labels (list of str): A list containing the labels for each dataset.
        """

        # If labels are not provided, use default labels
        if labels is None:
            labels = [f"Dataset {i+1}" for i in range(len(data_list))]

        # Determine global x and y limits
        all_x = [x for dataset in data_list for x in dataset[0]]
        all_y = [y for dataset in data_list for y in dataset[1]]

        x_min, x_max = min(all_x), max(all_x)
        y_min, y_max = min(all_y), max(all_y)

        # Calculate the number of rows and columns for subplots
        n = len(data_list)
        ncols = int(np.ceil(np.sqrt(n)))
        nrows = int(np.ceil(n / ncols))

        fig, axes = plt.subplots(nrows, ncols, figsize=(3*ncols, 3*nrows))

        # Make sure axes is always a 2D array, even when n is 1
        if n == 1:
            axes = np.array([axes])

        for ax, (x, y), label, annotation in zip(axes.ravel(), data_list, labels, annotations):
            ax.scatter(x, y, label='')
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(y_min, y_max)
            #ax.legend(fontsize=10)
            ax.set_title(label)

            ax.text(0.05, 0.1,
                    annotation,
                    transform=ax.transAxes,
                    verticalalignment='top')
            # Annotate each data point

        # Turn off any remaining unused subplots
        for i in range(n, nrows*ncols):
            axes.ravel()[i].axis('off')

        plt.tight_layout()
        plt.show()
