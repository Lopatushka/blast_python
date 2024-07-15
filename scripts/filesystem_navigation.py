import os
import shutil

def create_dir(new_directory):
    if not os.path.exists(new_directory):
        os.makedirs(new_directory)

def move_files(dir, files):
    for file in files:
        try:
            if os.path.isfile(file):
                destination_file = os.path.join(dir, os.path.basename(file))
                if os.path.exists(destination_file):
                    os.remove(destination_file)
                    print(f"Existing file {destination_file} removed.")
                shutil.move(file, dir)
        except shutil.Error as e:
            print(f"Error moving file {file} to {dir}: {e}")
        except OSError as e:
            print(f"OS error occurred while moving file {file}: {e}")
