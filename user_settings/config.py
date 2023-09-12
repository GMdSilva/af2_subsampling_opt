""" User-configured values """
from src.utilities.utilities import is_jupyter
import ipywidgets as widgets
from IPython.display import display, clear_output


SYSTEM_NAME = 'abl_wt'
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

def build_configs_dict():
    global SYSTEM_NAME,\
        PREDICTION_ROOT,\
        FIRST_RESIDUE,\
        REINDEX,\
        BUILD_MSA,\
        RUN_AF2,\
        OPTIMIZE_PARAMETERS,\
        TEST_MUTANTS,\
        JACKHMMER_PATH,\
        COLABFOLDBATCH_PATH,\
        MSA_PATHS

    configs_dict = {
        "SYSTEM_NAME": SYSTEM_NAME,
        "IS_JUPYTER": IS_JUPYTER,
        "PREDICTION_ROOT": PREDICTION_ROOT,
        "REINDEX": REINDEX,
        "FIRST_RESIDUE":  FIRST_RESIDUE,
        "MSA_PATHS": MSA_PATHS,
        "BUILD_MSA": BUILD_MSA,
        "JACKHMMER_PATH": JACKHMMER_PATH,
        "RUN_AF2": RUN_AF2,
        "COLABFOLDBATCH_PATH": COLABFOLDBATCH_PATH,
        "OPTIMIZE_PARAMETERS": OPTIMIZE_PARAMETERS,
        "TEST_MUTANTS": TEST_MUTANTS
    }
    return configs_dict

def update_config(configs_dict):
    global SYSTEM_NAME,\
        PREDICTION_ROOT,\
        FIRST_RESIDUE,\
        REINDEX,\
        BUILD_MSA,\
        RUN_AF2,\
        OPTIMIZE_PARAMETERS,\
        TEST_MUTANTS,\
        JACKHMMER_PATH,\
        COLABFOLDBATCH_PATH,\
        MSA_PATHS

    SYSTEM_NAME = configs_dict['SYSTEM_NAME']
    PREDICTION_ROOT = configs_dict['PREDICTION_ROOT']
    FIRST_RESIDUE = configs_dict['FIRST_RESIDUE']
    REINDEX = configs_dict['REINDEX']
    BUILD_MSA = configs_dict['BUILD_MSA']
    RUN_AF2 = configs_dict['RUN_AF2']
    OPTIMIZE_PARAMETERS = configs_dict['OPTIMIZE_PARAMETERS']
    TEST_MUTANTS = configs_dict['TEST_MUTANTS']
    JACKHMMER_PATH =configs_dict['JACKHMMER_PATH']
    COLABFOLDBATCH_PATH = configs_dict['COLABFOLDBATCH_PATH']
    MSA_PATHS = configs_dict['MSA_PATHS']



