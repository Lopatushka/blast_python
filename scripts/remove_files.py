import os
import subprocess
import glob
import warnings
from .get_files import list_of_files

def remove_files_by_pattern(directory, pattern):
    # List all files in the directory
    fasta_files = list_of_files(directory, "fa")
    matching_files = [file for file in fasta_files if pattern in file]
    if matching_files:
        for file in matching_files:
            subprocess.run(['rm', file], check=True)
    else:
         warnings.warn(f"There is no files in the directory {directory} with pattern {pattern} to delete.")

def remove_file_by_pattern_old(dir, pattern):
    path = os.path.join(dir, pattern)
    matching_files = glob.glob(path)
    if matching_files:
        for file in matching_files:
            subprocess.run(['rm', file], check=True)
    else:
        pass

def remove_file(file_path):
    if os.path.exists(file_path):
        try:
            subprocess.run(['rm', file_path], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error occurred: {e}")
    else:
        pass

# dir = "../sanger_seq/240822"
# pattern = ".fq"
# path = os.path.join(dir, pattern)
# matching_files = glob.glob(path)
# print(matching_files)