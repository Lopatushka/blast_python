import os
import subprocess
import glob
import warnings
from .get_files import list_of_files

def _remove_files_by_pattern(directory, extention, pattern):
    """
    EXTRA FUNCTION
    """
    fasta_files = list_of_files(directory, extention)
    matching_files = [file for file in fasta_files if pattern in file]
    if matching_files:
        for file in matching_files:
            subprocess.run(['rm', file], check=True)
    else:
         warnings.warn(f"There is no files in the directory {directory} with pattern {pattern} to delete.")

def remove_file(file_path):
    if os.path.exists(file_path):
        try:
            subprocess.run(['rm', file_path], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error occurred: {e}")
    else:
        pass

def remove_files_with_extension(dir, extension):
    pattern  = os.path.join(dir, f'*.{extension}')
    matching_files = glob.glob(pattern)
    if matching_files:
        for file_path  in matching_files:
            subprocess.run(['rm', file_path], check=True)
    else:
        print(f"There is no files with the extension .{extension} in the directory {dir}")