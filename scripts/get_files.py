import os
import argparse
import warnings

def is_valid_directory(path):
    """
    Check if the directory exists."""
    if not os.path.isdir(path):
        raise argparse.ArgumentTypeError(f"Directory '{path}' does not exist.")
    return path

def check_range(min_value, max_value):
    """Return a function that checks if an integer is within the desired range."""
    def validate(value):
        try:
            value = int(value)
        except ValueError as e:
            raise argparse.ArgumentTypeError(e)
        if value < min_value or value > max_value:
            raise argparse.ArgumentTypeError(f"Value '{value}' is outside the range [{min_value}, {max_value}].")
        return value
    return validate

def check_blast_database(database):
    try:
        # Check if dir exists
        dirname = os.path.dirname(database)
        if not os.path.isdir(dirname):
            raise argparse.ArgumentTypeError(f"Path '{database}' to blastn database is not found")
        
        # Check taxonomy_files are located in blast database dir
        list_of_files = os.listdir(dirname)
        taxonomy_files = ['taxdb.btd', 'taxdb.bti', 'taxonomy4blast.sqlite3']
        for taxonomy_file in taxonomy_files:
            if taxonomy_file not in list_of_files:
                raise argparse.ArgumentTypeError(f"{taxonomy_file} is not found in the blastn database directory '{dirname}'")
            
        # Check validity of database name     
        database_name = database.split("/")[-1]
        for taxonomy_file in taxonomy_files:
            list_of_files.remove(taxonomy_file)
        if sum([file[:-4] == database_name for file in list_of_files]) == 0:
            raise argparse.ArgumentTypeError(f"Wrong blastn database name: '{database_name}'")
        return database
    except argparse.ArgumentTypeError:
        raise # Re-raise exceptions above 
    except Exception as e:
        raise argparse.ArgumentTypeError(f"Unexpected error occured during --database argument checking: {e}")

def parse_arguments():
    parser = argparse.ArgumentParser(description="Reads trimming, build consensus and perform blastn.")
    parser.add_argument('--directory', '-d', type=is_valid_directory, required=True, help='Directory containing .ab1 files')
    parser.add_argument('--parsing_mode', '-pm', type=str, default="auto", choices=["auto", "manual"], help='Filenames parsing mode. Options: auto (default), manual')
    parser.add_argument('--parsing_patterns', '-pp', nargs='*', type=str, default=list(), help='Unique patterns in filenames of .ab1 files for manual consensus building')
    parser.add_argument('--blastn_mode', '-bm', type=str, default="auto", choices=["auto", "manual"], help='Balstn search mode. Options: auto (default), manual')
    parser.add_argument('--consensus_patterns', '-cp', nargs='*', type=str, default=list(), help='Unique patterns in filenames of consensuses for manual blastn search')
    
    parser.add_argument('--consensus_quality', '-cq', type=check_range(1, 100), default=10,
                        help="Maximum percetage of 'N' or degenerate nucleotides in consensus. Permissiable value range is [1, 100]. Default value is 10.")
    parser.add_argument('--nthreads', '-nt', type=check_range(1, 100), default=4,
                        help='Number of threads for blastn search. Permissiable value range is [1, 100]. Default value is 4.')
    parser.add_argument('--trimming_quality', '-tq', type=check_range(1, 100), default=15,
                        help='Trimming quality. Permissiable value range is [1, 100]. Default value is 15.')
    parser.add_argument('--minlength', '-ml', type=check_range(50, 1000), default=50,
                        help='Minimum length of consensus. Permissiable value range is [1, 1000]. Default value is 50.')
    parser.add_argument('--database', '-db', type=check_blast_database, required=True, help='Blastn database')
    parser.add_argument('--qcovus', '-qc', type=check_range(1, 100), default=80,
                        help='Query coverage threshold. Permissiable value range is [1, 100]. Default value is 80.')
    parser.add_argument('--pident', '-pi', type=check_range(1, 100), default=95,
                        help='Percent of identity threshold. Permissiable value range is [1, 100]. Default value is 95.')
    return parser.parse_args()

def list_of_files(dir, extension):
        '''
        Description:
            This function returns a list of full file paths in a given directory that have a specific file extension.
            The function iterates over all files in the directory and filters out those that do not match the specified extension.

        Parameters:
            dir (str): The path to the directory where the search will be performed.
            extension (str): The file extension to filter files (e.g., 'txt', 'jpg'). Files with this extension will be returned.

        Returns:
            list: A list containing full file paths of the files with the specified extension found in the directory.
                If the directory is not found or other errors occur, an empty list is returned.
        '''
        try:
            files = []
            for item in os.listdir(dir):
                path=os.path.join(dir, item)
                if os.path.isfile(path):
                    files.append(path)
            subset_of_files=[file for file in files if f'.{extension}' in file] # full paths to .ext files
            return subset_of_files
        except FileNotFoundError as e:
            warnings.warn(f"Directory '{dir}' is not found: {e}")
            return []
        except PermissionError as e:
            warnings.warn(f"Permission denied for directory '{dir}': {e}")
            return []
        except Exception as e:
            warnings.warn(f"An unexpected error occurred while listing files in directory '{dir}': {e}")
            return []

def list_of_files_by_pattern(dir, extension, patterns):
    '''
    Description:
        This function returns a list of files from a specified directory that match both a given file extension and 
        a set of patterns. It builds on top of the `list_of_files` function, which retrieves all files with the specified 
        extension. Then, it filters these files to return only those whose filenames contain one or more of the given patterns.

    Parameters:
        dir (str): The path to the directory where the search will be performed.
        extension (str): The file extension to filter files (e.g., 'txt', 'csv'). Only files with this extension are considered.
        patterns (list): A list of strings (patterns) to match against filenames. Files that contain any of these patterns 
                         in their filename will be included in the returned list.

    Returns:
        list: A list of file paths that match both the specified extension and the provided patterns.
              If no files match, an empty list is returned.
    '''
    files = list_of_files(dir, extension=extension) # get all files with reqired extension
    files_to_process = []
    for pattern in patterns:
        for file in files:
            if pattern in file:
                files_to_process.append(file)
    return files_to_process


def filename_parsing(file):
    try:
        sample_types = ["C", "CS", "F", "A", "P"]
        #primer_names = []
        
        # Get filename without extension
        filename = file.split("/")[-1][:-4]

        # Split filename by _ and check its length
        filename_splitted = filename.split("_")

        # Check filename splitting
        # Length array after splitting
        if len(filename_splitted) < 4:
            raise ValueError(f"Error parsing filename {file}. "
                                f"Filename has to be started with Plate-YY-MM-DD and contain the following information: "
                                f"sample_type, sample_name, primer_name and primer orientation. "
                                f"For detailes see READ.me")

        # Plate-YY-MM-DD
        if (len(filename_splitted[0].split("-")) != 4) | (filename_splitted[0].split("-")[0] != "Plate"):
            raise ValueError(f"Error parsing filename {file}: filename has to be started with 'Plate-YY-MM-DD'."
                            f"For detailes see READ.me")
        
        # Get sample type and check it
        sample_type = filename_splitted[1]
        if sample_type not in sample_types:
            raise ValueError(f"Error parsing filename {file}: wrong sample_type: {sample_type}. "
                            f"Possible options for sample_type: {sample_types}. For detailes see READ.me")
    
        # Get sample name
        sample_name = filename_splitted[2]

        # Get primer name
        primer = filename_splitted[3]
        primer_splitted = primer.split("-")
        primer_name = primer_splitted[0]
        primer_orientation = primer_splitted[1]
        sample_name_primer = sample_name + "_" + primer_name

        # Check primer name
        if len(primer_splitted) != 2:
            raise ValueError(f"Error parsing filename {file}: wrong primer name: {primer}. "
                            f"For detailes see READ.me")
        elif ("F" not in primer_splitted[1]) & ("f" not in primer_splitted[1]) & ("R" not in primer_splitted[1]) & ("r" not in primer_splitted[1]):
            raise ValueError(f"Error parsing filename {file}: wrong primer name: {primer}. "
                            f"You need to specify primer's orientartion. "
                            f"For detailes see READ.me")
        #elif primer_splitted[0] not in primer_names:
            #raise ValueError(f"Error parsing filename {file}: wrong primer name {primer}.\
                                # For detailes see READ.me")
        
        # Get path to dir
        dir = "/".join(file.split("/")[:-1])
        # Check dir and parsing prosedure in general
        if not os.path.isdir(dir):
            raise ValueError(f"Error parsing filename {file}: directory {dir} doesn't exist")
        
        return {'filename': filename,
                    'sample_type': sample_type,
                    'sample_name': sample_name,
                    'primer': primer,
                    'sample_name_primer':sample_name_primer,
                    'path': file,
                    'dir': dir}
    except ValueError as e:
        print(f"Error occured: {e}")