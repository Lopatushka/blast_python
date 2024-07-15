import os
from pathlib import Path

def check_1_file_exists(file_path):
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    
def check_files_exist(file_paths_list):
    for filepath in file_paths_list:
        check_1_file_exists(filepath)

def is_empty(data):
    if not data:
        raise ValueError(f"Variable {data} is empty. Exiting the program.")
    
def not_empty_file(file_path):
    return Path(file_path).stat().st_size != 0