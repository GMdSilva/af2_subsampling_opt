""" User-configured values """

from user_settings.new_config import load_config

CONFIG = load_config()
PREFIX = CONFIG['PREFIX']
REINDEX = CONFIG['REINDEX']
FIRST_RESIDUE = CONFIG['FIRST_RESIDUE']
