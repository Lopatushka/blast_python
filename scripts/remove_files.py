import os
import subprocess
import glob

def remove_file_by_pattern(dir, pattern):
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