import os
import warnings
class FileOverwrittenWarning(UserWarning):
    pass

def check_1_file_exists(file_path):
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    
def check_files_exist(file_paths_list):
    for filepath in file_paths_list:
        check_1_file_exists(filepath)
    
def file_already_exists(file_path):
    try:
        if os.path.isfile(file_path):
            if os.path.getsize(file_path) != 0:
                warnings.warn(f"File {file_path} will be overwritten", FileOverwrittenWarning)
    except Exception as e:
        print(f"There is an error: {e} during file checking {file_path}")

def is_file_exists(path):
    """Check if the file exists and not empty."""
    if (os.path.isfile(path) and os.path.getsize(path) != 0):
        return True
    return False