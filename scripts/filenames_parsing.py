import argparse
import os
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser(description="Reads trimming, build consensus and perform blastn.")
    parser.add_argument('--input_directory', type=str, default=".", help='Directory containing .ab1 files')
    parser.add_argument('--output_directory', type=str, default=".", help='Directory to save results')
    parser.add_argument('--files_processing_mode', type=str, default="auto", help='Files processing mode. Options: auto (default), manual')
    parser.add_argument('--filenames_patterns', type=list, default=list())
    parser.add_argument('--blastn_mode', type=str, default="auto", help='Balstn search mode. Options: auto (default), manual')
    parser.add_argument('--consensus_patterns', type=list, default=list())
    parser.add_argument('--nthreads', type=int, default=4, help='Number of threads for blastn search')
    parser.add_argument('--trimming_quality', type=int, default=15, help='Trimming quality')
    parser.add_argument('--minlength', type=int, default=50, help='Minimum length of consensus')
    parser.add_argument('--database', type=str, default=".", help='Blastn database')
    parser.add_argument('--qcovus', type=int, default=80, help='Query coverage threshold')
    parser.add_argument('--pident', type=int, default=95, help='Percent of identity threshold')
    return parser.parse_args()

def filenames_parsing():
    args = parse_arguments()
    dir = args.input_directory
    dir = "../blast/data/101424"

    files = []
    for item in os.listdir(dir):
        path=os.path.join(dir, item)
        if os.path.isfile(path):
            files.append(path)
    ab1_files=[file for file in files if ".ab1" in file]
    filenames = [file.split("/")[-1][:-4] for file in ab1_files] # without extension
    sample_types = [file.split("_")[1] for file in ab1_files]
    sample_names = [file.split("_")[2] for file in ab1_files]
    primers = [file.split("_")[3] for file in ab1_files]

    return zip(filenames, sample_types, sample_names, primers)
    
def main():
    for f,t,n,p in filenames_parsing():
        print(f,t,n,p)

if __name__ == "__main__":
    main()