""" User-configured values """
from angra.utilities.utilities import is_jupyter
import ipywidgets as widgets
from IPython.display import display, clear_output


SYSTEM_NAME = 'abl_wt'
IS_JUPYTER = is_jupyter()
PREDICTION_ROOT = "."
REINDEX = False,
FIRST_RESIDUE = 0
MSA_PATHS = ''
BUILD_MSA = False
JACKHMMER_PATH = ''
RUN_AF2 = False
COLABFOLDBATCH_PATH = ''
OPTIMIZE_PARAMETERS = True
TEST_MUTANTS = True


# Create widgets for each of the parameters
if is_jupyter():
        system_name_widget = widgets.Text(value='',
                                          description='Label for predictions/results')
        prediction_root_widget = widgets.Text(value='.',
                                              description='Path to save predictions/results to')
        reindex_residues_widget = widgets.Checkbox(value=False,
                                                   description='Reindex Residues (default is 1 indexed)')
        first_residue_widget = widgets.Text(value='1',
                                            description='If reindexing, first residue index')
        msa_paths_widget = widgets.Text(value=f'results/msas/',
                                        description='Path to saved MSAs')
        build_msa_widget = widgets.Checkbox(value=False,
                                            description='Use custom MSAs')
        jackhmmer_path_widget = widgets.Text(value='tools/jackhmmer/',
                                             description='Path to jackhmmer')
        run_af2_widget = widgets.Checkbox(value=False,
                                          description='Make predictions using AF2')
        colabfoldbatch_path_widget = widgets.Text(value='tools/colabfoldbatch/',
                                                  description='Path to colabfoldbatch')
        optimize_parameters_widget = widgets.Checkbox(value=True,
                                                      description='Optimize AF2 predictions')
        test_mutants_widget = widgets.Checkbox(value=True,
                                               description=f'Test mutants of')


        # Display the widgets

        display(system_name_widget, prediction_root_widget, reindex_residues_widget, first_residue_widget,
                msa_paths_widget, build_msa_widget, jackhmmer_path_widget, run_af2_widget, colabfoldbatch_path_widget,
                optimize_parameters_widget, test_mutants_widget)


        def update_config():
            global SYSTEM_NAME, PREDICTION_ROOT, FIRST_RESIDUE, REINDEX, BUILD_MSA, RUN_AF2, OPTIMIZE_PARAMETERS, TEST_MUTANTS, JACKHMMER_PATH, COLABFOLDBATCH_PATH, MSA_PATHS

            SYSTEM_NAME = system_name_widget.value
            PREDICTION_ROOT = prediction_root_widget.value
            FIRST_RESIDUE = first_residue_widget.value
            REINDEX = reindex_residues_widget.value
            BUILD_MSA = build_msa_widget.value
            RUN_AF2 = run_af2_widget.value
            OPTIMIZE_PARAMETERS = optimize_parameters_widget.value
            TEST_MUTANTS = test_mutants_widget.value
            JACKHMMER_PATH = jackhmmer_path_widget.value
            COLABFOLDBATCH_PATH = colabfoldbatch_path_widget.value
            MSA_PATHS = msa_paths_widget.value


        update_config_button = widgets.Button(description="Update Config")
        update_config_button.on_click(lambda b: update_config())
        display(update_config_button)

        # Display the widgets

        display(system_name_widget, prediction_root_widget, reindex_residues_widget, first_residue_widget,
                msa_paths_widget, build_msa_widget, jackhmmer_path_widget, run_af2_widget, colabfoldbatch_path_widget,
                optimize_parameters_widget, test_mutants_widget)


