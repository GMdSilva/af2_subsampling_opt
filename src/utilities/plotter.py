"""
Defines a class for plotting
"""

from typing import Dict, List, Union, Tuple
import math
import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

from user_settings.config import PREDICTION_ROOT, IS_JUPYTER, SYSTEM_NAME


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
        plt.close()
        self.data = data
        self.prefix = prefix

        plt.rcParams['xtick.labelsize'] = 16  # x ticks
        plt.rcParams['ytick.labelsize'] = 16  # y ticks
        plt.rcParams['axes.labelsize'] = 18  # x and y labels
        plt.rcParams['axes.titlesize'] = 20  # title of each subplot
        plt.rcParams['figure.titlesize'] = 20  # overall figure title
        plt.rcParams['legend.fontsize'] = 16  # figure legend fontsize

    def plot_rmsf_peaks(self,
                        subsampling_results: List[Dict],
                        peak_ranges: List[List[Tuple[int, int]]],
                        save_plot: bool = True) -> None:
        """
        Plots peaks for the given results.

        Parameters:
        - subsampling_results: A list of data dictionaries.
        - peak_ranges: A list of peak range lists corresponding to each data set.
        - save_plot (bool): If plot is to be saved to disk or not.
        """
        plt.close()

        colors = ['red', 'blue', 'green', 'yellow', 'purple',
                  'orange', 'teal', 'fuchsia', 'pink']

        num_data_sets = len(subsampling_results)

        # Calculate the optimal number of rows and columns for the subplots.
        num_cols = 2
        num_rows = math.ceil(num_data_sets / num_cols)

        _, axs = plt.subplots(num_rows, num_cols, figsize=(7 * num_cols, 3 * num_rows),
                              sharex='col')
        if num_data_sets == 1:
            axs = [axs]
        elif num_data_sets <= num_cols:
            axs = [ax for ax in axs]
        else:
            axs = [ax for sublist in axs for ax in sublist]  # Flatten the axs list.

        for ax in axs[num_data_sets:]:  # Turn off the axes for any unused subplots.
            ax.axis('off')

        for idx, (data, ax) in enumerate(zip(subsampling_results, axs)):
            ax.plot(data['residues'], data['results'],
                    marker=None,
                    label=data['trial'],
                    color='grey')

            # Only label the y-axis for the leftmost plots.
            if idx % num_cols == 0:
                ax.set_ylabel('C-Alpha RMSF (A)',
                              weight='bold')

            ax.set_xlabel('Residue #',
                          weight='bold')
            ax.legend()
            for i, (start, end) in enumerate(peak_ranges):
                color = colors[i % len(colors)]
                ax.plot(data['residues'][start:end],
                        data['results'][start:end],
                        marker='*',
                        color=color,
                        linestyle='-',
                        markersize=plt.rcParams['lines.markersize'] / 2)
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.15,
                            wspace=0.1)  # Adjust this value for desired vertical spacing.
        save_path = os.path.join(PREDICTION_ROOT,
                                 'results',
                                 'plots',
                                 f"{self.prefix}_rmsf_peak_detection.png")
        if save_plot:
            plt.savefig(save_path)
        if IS_JUPYTER:
            plt.show()

    def plot_kde_with_modes(self,
                            residue_ranges: List[Tuple[int, int]],
                            save_plot=True) -> None:
        """
        Plots Kernel Density Estimation (KDE) with mode peaks for a list of data sets.

        This function generates multiple KDE plots based on the provided data,
        highlights the mode peaks in red, and saves the combined plot to the specified path.

        Parameters:
        - data_sets (list): A list of data dictionaries.
        - peak_ranges (list): A list of tuples each representing
                a range of peaks for corresponding data set.

        Returns:
        None. The function saves the generated multiplot to the specified path.
        """
        plt.close()
        num_data_sets = len(self.data)
        num_cols = 2
        num_rows = math.ceil(num_data_sets / num_cols)

        _, axs = plt.subplots(num_rows, num_cols, figsize=(7 * num_cols, 3 * num_rows))
        if num_data_sets == 1:
            axs = [axs]
        elif num_data_sets <= num_cols:
            axs = [ax for ax in axs]
        else:
            axs = [ax for sublist in axs for ax in sublist]  # Flatten the axs list.

        for i, (data, ax) in enumerate(zip(self.data, axs)):
            ax.plot(data["x_values"],
                    data["density"],
                    label=data["trial"])
            ax.scatter(data["x_values"][data["peaks"]],
                       data["density"][data["peaks"]],
                       color='red')
            ax.axhline(y=data['threshold'],
                       color='green',
                       linestyle='--')
            ax.legend()

            ax.set_xlabel('Residue #',
                          weight='bold')

            # Only label the y-axis for the left-most subplots
            if i % num_cols != 0:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel('Density', weight='bold')

            ax.set_title(f'Residue Range {str(residue_ranges)}')

        # Hide any unused subplots.
        for i in range(num_data_sets, num_rows * num_cols):
            axs[i].axis('off')

        plt.tight_layout()

        if IS_JUPYTER:
            plt.show()

        if save_plot:
            save_path = os.path.join(PREDICTION_ROOT,
                                     'results',
                                     'plots',
                                     f"{self.prefix}_{residue_ranges[0]}_{residue_ranges[1]}"
                                     f"_kde_modes.png")
            plt.savefig(save_path)

    def plot_mut_analysis(self,
                          x_values: Dict,
                          mut_list: Dict,
                          known_ground_effects: bool = False,
                          known_alt1_effects: bool = False) -> None:
        """
        Generate a bar plot based on mutation analysis.

        Parameters:
        - x_values (dict): A dictionary containing the key-value pairs for the x-axis.
        - mut_list (dict): A dictionary containing mutation details.
        - known_ground_effects (bool): If effects on ground state are known or not.
        - known_alt1_effects (bool): If effects on alt1 state are known or not.

        Returns:
        None. A bar plot is saved in the designated path.
        """
        plt.close()
        ground_pop_diffs = [d[self.prefix][list(x_values.keys())[0]]
                            for d in self.data for self.prefix in d]
        labels = [entry['label'] for entry in mut_list.values()]
        effects_ground = [entry['effect']['ground_pop']
                          for entry in mut_list.values()]
        effects_alt1 = [entry['effect']['alt1_pop']
                        for entry in mut_list.values()]

        colors = ['red' if ground_pop_diff < 0
                  else ('grey' if ground_pop_diff == 0 else 'blue')
                  for ground_pop_diff in ground_pop_diffs]

        if known_ground_effects:
            if x_values == 'ground_pop_diff' or 'ground_pop_test':
                colors = ['red' if effect == '-'
                          else ('grey' if effect == "ref"
                                else 'blue') for effect in effects_ground]
        if known_alt1_effects:
            if x_values == 'alt1_pop_diff' or 'alt1_pop_test':
                colors = ['red' if effect == '-'
                          else ('grey' if effect == "ref"
                                else 'blue') for effect in effects_alt1]
        # Plotting
        _, ax = plt.subplots()
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
                                 f"{SYSTEM_NAME}_{list(x_values.keys())[0]}.png")
        plt.savefig(save_path)

    def plot_multiple_scatter(self,
                              data_list,
                              save_path=None,
                              labels=None,
                              annotations=None,
                              colors=None,
                              colorbar_label=None,
                              fontsize=10):
        """
        Creates a separate scatter plot for each dataset,
            organized in an automatically determined grid.

        Args:
        - data_list (list of tuples): A list where each element
                is a tuple containing two lists (x_values, y_values).
        - labels (list of str): A list containing the labels for each dataset.
        - annotations (list of str): A list containing the annotations for each dataset.
        - colors (list of list): A list containing color values
            for each data point in the datasets.
        - colorbar_label (str): Label for the colorbar.
        """
        plt.close()
        # If labels are not provided, use default labels
        if labels is None:
            labels = [f"Dataset {i + 1}" for i in range(len(data_list))]

        # Determine global x and y limits
        all_x = [x for dataset in data_list for x in dataset[0]]
        all_y = [y for dataset in data_list for y in dataset[1]]

        x_min, x_max = min(all_x) - .5, max(all_x) + .5
        y_min, y_max = min(all_y) - .5, max(all_y) + .5

        # Calculate the number of rows and columns for subplots
        n = len(data_list)
        ncols = int(np.ceil(np.sqrt(n)))
        nrows = int(np.ceil(n / ncols))

        # Use gridspec to organize the layout
        fig = plt.figure(figsize=(3 * ncols + 2, 3 * nrows))
        gs = gridspec.GridSpec(nrows,
                               ncols + 1,
                               width_ratios=[1] * ncols + [0.1],
                               wspace=0.25,
                               hspace=0.35)

        # Use a common color map
        cmap = plt.get_cmap('viridis')

        # If using global color limits, determine them here
        color_min = min([min(c) for c in colors])
        color_max = max([max(c) for c in colors])

        sc = None  # To hold the scatter plot object for the colorbar reference

        for i, ((x, y), label, annotation, color) in enumerate(zip(data_list,
                                                                   labels,
                                                                   annotations,
                                                                   colors)):
            ax = fig.add_subplot(gs[i // ncols, i % ncols])
            sc = ax.scatter(x,
                            y,
                            c=color,
                            cmap=cmap,
                            vmin=color_min,
                            vmax=color_max,
                            s=50,
                            edgecolor='black',
                            linewidth=1)
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(y_min, y_max)
            ax.set_title(label, fontsize=15, weight='bold')
            ax.text(0.01,
                    .1,
                    annotation,
                    transform=ax.transAxes,
                    verticalalignment='top',
                    weight='bold',
                    fontsize=fontsize)

        # Add a single colorbar outside the grid using the space defined by GridSpec
        if colorbar_label:
            cbar_ax = fig.add_subplot(gs[:, -1])
            cbar = fig.colorbar(sc, cax=cbar_ax, orientation='vertical')

        # Set colorbar label if provided
            cbar.set_label(colorbar_label, fontsize=15)
        if save_path:
            save_path = os.path.join(PREDICTION_ROOT,
                                     'results',
                                     'plots',
                                     f"{self.prefix}_mutants_ground_vs_alt1.png")
            plt.savefig(save_path)
        if IS_JUPYTER:
            plt.show()

    def plot_bar_for_state(self, state, mut_list, known_ground_effects = False):
        plt.close()
        labels = list(self.data.keys())
        values = []
        for label in labels:
            values.append(self.data[label][state])

        colors = ['red' if value < 0
                  else ('grey' if value == 0 else 'blue')
                  for value in values]

        if known_ground_effects:
            if state == 'Ground':
                effects_ground = [entry['effect']['ground_pop']
                                  for entry in mut_list.values()]

                colors = ['red' if effect == '-'
                          else ('grey' if effect == "ref"
                                else 'blue') for effect in effects_ground]

        # Plotting
        _, ax = plt.subplots()
        ax.barh(labels,
                values,
                color=colors)
        ax.set_xlabel(f"{state} State Pop. Δ (%)")
        ax.set_ylabel('')
        ax.set_title('')

        plt.tight_layout()
        save_path = os.path.join(PREDICTION_ROOT,
                                 'results',
                                 'plots',
                                 f"{SYSTEM_NAME}_{state}_diff_2danalysis.png")
        plt.savefig(save_path)