import os
from pathlib import Path
import warnings
from get_files import filename_parsing

def check_1_file_exists(file_path):
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    
def check_files_exist(file_paths_list):
    for filepath in file_paths_list:
        check_1_file_exists(filepath)

def is_empty_variable(var, var_name):
    if not var:
        raise ValueError(f"Variable {var_name} is empty. Exiting the program.")
    
def not_empty_file(file_path):
    return Path(file_path).stat().st_size != 0

class ArgumentError(Exception):
    """Custom exception for invalid arguments."""
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

def check_blast_database(database):
    dir = os.path.dirname(database)
    if not os.path.isdir(dir):
        raise IOError(f"Path {dir} to blastn database is not found")
    
    list_of_files = os.listdir(dir)

    # Check taxonomy_files are located in blast database dir
    taxonomy_files = ['taxdb.btd', 'taxdb.bti', 'taxonomy4blast.sqlite3']
    for taxonomy_file in taxonomy_files:
        if taxonomy_file not in list_of_files:
            raise IOError(f"{taxonomy_file} is not found in blastn database {database}")
    database_name = database.split("/")[-1]
    
    for taxonomy_file in taxonomy_files:
        list_of_files.remove(taxonomy_file)
    
    if sum([file[:-4] == database_name for file in list_of_files]) == 0:
        raise ValueError(f"Wrong value of blastn database name: {database_name}")
    
    def check_ab1_filenames(directory):
        try:
            filenames = list_of_files(dir=directory, extension="ab1")
            for file in filenames:
                filename_parsing(file)
        except FileNotFoundError as e:
                warnings.warn(f"Directory '{dir}' is not found: {e}")
                return []
        except PermissionError as e:
                warnings.warn(f"Permission denied for directory '{dir}': {e}")
                return []
        except ValueError as e:
                warnings.warn(f"Error in file processing: {e}")         
        except Exception as e:
                warnings.warn(f"An unexpected error occurred: {e}")
                return []