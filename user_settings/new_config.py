""""Functions for editing and saving config file"""

import json

def set_config(prefix: str,
               path: str,
               reindex: bool = False,
               first_residue: int = 0) -> None:

    config_data = {
        "PREFIX": prefix,
        "PATH": path,
        "REINDEX": reindex,
        "FIRST_RESIDUE": first_residue,
    }

    with open("config.json", "w", encoding="utf-8") as config_file:
        json.dump(config_data, config_file, indent=4)


def load_config(filename="config.json"):
    try:
        with open(filename, "r", encoding="utf-8") as config_file:
            config_data = json.load(config_file)
            return config_data
    except FileNotFoundError:
        print(f"{filename} not found.")
        return {}
    except json.JSONDecodeError:
        print(f"Error decoding {filename}. Please ensure it's valid JSON.")
        return {}
