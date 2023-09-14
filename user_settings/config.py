""" User-configured values """
import ipywidgets as widgets
from IPython.display import display, clear_output

def is_jupyter():
    try:
        # The `get_ipython` function is available in Jupyter environments, including IPython.
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type
    except NameError:
        return False      # Probably standard Python interpreter

SYSTEM_NAME = 'gmcsf_wt'
IS_JUPYTER = is_jupyter()
PREDICTION_ROOT = "."
REINDEX = False
FIRST_RESIDUE = 1
MSA_PATHS = ''
BUILD_MSA = False
JACKHMMER_PATH = ''
RUN_AF2 = False
COLABFOLDBATCH_PATH = ''
OPTIMIZE_PARAMETERS = True
TEST_MUTANTS = True

# MUTANT_DATA = {
#     "abl_wt": {"label": "Abl1", "effect": {"ground_pop":"ref", "alt1_pop":"ref"}, "rank": 4},
#     "abl_E255KandT315I": {"label": "Abl1 E255K+T315I", "effect": {"ground_pop":"+", "alt1_pop":"-"}, "rank": 3},
#     "abl_E255V": {"label": "Abl1 E255V", "effect": {"ground_pop":"+", "alt1_pop":"-"}, "rank": 0},
#     "abl_E255VandT315I": {"label": "Abl1 E255V+T315I", "effect": {"ground_pop":"+", "alt1_pop":"-"}, "rank": 2},
#     "abl_F382L": {"label": "Abl1 F382L", "effect": {"ground_pop":"+", "alt1_pop":"-"}, "rank": 1},
#     "abl_F382V": {"label": "Abl1 F382V", "effect": {"ground_pop":"-", "alt1_pop":"+"}, "rank": 5},
#     "abl_F382Y": {"label": "Abl1 F382Y", "effect": {"ground_pop":"-", "alt1_pop":"+"}, "rank": 6},
#     "abl_I2M": {"label": "Abl1 I2M", "effect": {"ground_pop":"-", "alt1_pop":"+"}, "rank": 7},
#     "abl_L301I": {"label": "Abl1 L301I", "effect": {"ground_pop":"-", "alt1_pop":"+"}, "rank": 8},
#     "abl_M290L": {"label": "Abl1 M290L", "effect": {"ground_pop":"-", "alt1_pop":"+"}, "rank": 9},
#     "abl_M290LandL301I": {"label": "Abl1 M290L+301I", "effect": {"ground_pop":"-", "alt1_pop":"+"}, "rank": 10}
# }

MUTANT_DATA = {
    "gmcsf_wt": {"label": "GMCSF", "effect": {"ground_pop":"ref", "alt1_pop":"ref"}, "rank": 0},
    "gmcsf_h15y": {"label": "GMCSF H15Y", "effect": {"ground_pop":"-", "alt1_pop":"+"}, "rank": 1},
    "gmcsf_h83y": {"label": "GMCSF H83Y", "effect": {"ground_pop":"-", "alt1_pop":"+"}, "rank": 2},
    "gmcsf_h15r": {"label": "GMCSF H15R", "effect": {"ground_pop":"-", "alt1_pop":"+"}, "rank": 3},
    "gmcsf_h87r": {"label": "GMCSF H87R", "effect": {"ground_pop":"-", "alt1_pop":"+"}, "rank": 4},
}

# MUTANT_DATA = {
#     "porb_wt": {"label": "PORB", "effect": {"ground_pop":"ref", "alt1_pop":"ref"}, "rank": 0},
#     "porb_d122n": {"label": "PORB D122N", "effect": {"ground_pop": "+", "alt1_pop": "-"}, "rank": 0},
#     "porb_e318n": {"label": "PORB E318N", "effect": {"ground_pop": "+", "alt1_pop": "-"}, "rank": 0},
#     "porb_s312e": {"label": "PORB S312E", "effect": {"ground_pop": "+", "alt1_pop": "-"}, "rank": 0}
# }


def build_configs_dict():
    global SYSTEM_NAME,\
        PREDICTION_ROOT,\
        FIRST_RESIDUE,\
        REINDEX,\
        JACKHMMER_PATH,\
        COLABFOLDBATCH_PATH,\
        MSA_PATHS,\
        MUTANT_DATA,\
        IS_JUPYTER

    configs_dict = {
        "SYSTEM_NAME": SYSTEM_NAME,
        "IS_JUPYTER": IS_JUPYTER,
        "PREDICTION_ROOT": PREDICTION_ROOT,
        "REINDEX": REINDEX,
        "FIRST_RESIDUE":  FIRST_RESIDUE,
        "MSA_PATHS": MSA_PATHS,
        "JACKHMMER_PATH": JACKHMMER_PATH,
        "RUN_AF2": RUN_AF2,
        "COLABFOLDBATCH_PATH": COLABFOLDBATCH_PATH,
        "MUTANT_DATA": MUTANT_DATA
    }
    return configs_dict

def update_config(configs_dict):
    global SYSTEM_NAME,\
        PREDICTION_ROOT,\
        FIRST_RESIDUE,\
        REINDEX,\
        JACKHMMER_PATH,\
        COLABFOLDBATCH_PATH,\
        MSA_PATHS,\
        MUTANT_DATA,\
        IS_JUPYTER

    SYSTEM_NAME = configs_dict['SYSTEM_NAME']
    PREDICTION_ROOT = configs_dict['PREDICTION_ROOT']
    FIRST_RESIDUE = configs_dict['FIRST_RESIDUE']
    REINDEX = configs_dict['REINDEX']
    JACKHMMER_PATH =configs_dict['JACKHMMER_PATH']
    COLABFOLDBATCH_PATH = configs_dict['COLABFOLDBATCH_PATH']
    MSA_PATHS = configs_dict['MSA_PATHS']
    MUTANT_DATA = configs_dict['MUTANT_DATA']
    IS_JUPYTER = configs_dict['IS_JUPYTER']


