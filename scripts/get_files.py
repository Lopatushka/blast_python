import os
import argparse
import warnings

def parse_arguments():
    parser = argparse.ArgumentParser(description="Reads trimming, build consensus and perform blastn.")
    parser.add_argument('--directory', '-d', type=str, help='Directory containing .ab1 files')
    parser.add_argument('--parsing_mode', '-pm', type=str, default="auto", choices=["auto", "manual"], help='Filenames parsing mode. Options: auto (default), manual')
    parser.add_argument('--parsing_patterns', '-pp', nargs='*', type=str, default=list(), help='Unique patterns in filenames of .ab1 files for manual consensus building')
    parser.add_argument('--blastn_mode', '-bm', type=str, default="auto", choices=["auto", "manual"], help='Balstn search mode. Options: auto (default), manual')
    parser.add_argument('--consensus_patterns', '-cp', nargs='*', type=str, default=list(), help='Unique patterns in filenames of consensuses for manual blastn search')
    parser.add_argument('--consensus_quality', '-cq', type=int, default=15, help="Maximum percetage of 'N' in consensus. Default value = 15")
    parser.add_argument('--nthreads', '-nt', type=int, default=4, help='Number of threads for blastn search')
    parser.add_argument('--trimming_quality', '-tq', type=int, default=15, help='Trimming quality. Default value = 15')
    parser.add_argument('--minlength', '-ml', type=int, default=50, help='Minimum length of consensus. Default value = 50')
    parser.add_argument('--database', '-db', type=str, default="/home/lopatushka/db/16S_ribosomal_RNA/16S_ribosomal_RNA", help='Blastn database')
    parser.add_argument('--qcovus', '-qc', type=int, default=80, help='Query coverage threshold. Default value = 80')
    parser.add_argument('--pident', '-pi', type=int, default=95, help='Percent of identity threshold. Default value = 95')
    return parser.parse_args()

def list_of_files(dir, extension):
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
    files = list_of_files(dir, extension=extension) # get all files with reqired extension
    files_to_process = []
    for pattern in patterns:
        for file in files:
            if pattern in file:
                files_to_process.append(file)
    return files_to_process

def filename_parsing(file):
    try:
        filename = file.split("/")[-1][:-4] # file name without extension
        sample_type = file.split("_")[1] # sample type: C - culture, P - plasmid etc.
        sample_name = file.split("_")[2] # sample name
        primer = file.split("_")[3] # primer
        dir = "/".join(file.split("/")[:-1]) # path to directory
        return {'filename': filename,
                'sample_type': sample_type,
                'sample_name': sample_name,
                'primer': primer,
                'path': file,
                'dir': dir}
    except IndexError as e:
        warnings.warn(f"Error parsing filename {file}: {e}")
        return {}
    except AttributeError as e:
        warnings.warn(f"Invalid input type for filename {file}: {e}")
        return {}
    except Exception as e:
        warnings.warn(f"An unexpected error occurred while parsing filename '{file}': {e}")
        return {}