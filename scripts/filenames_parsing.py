import argparse
import os
import subprocess
import numpy as np
from Bio import SeqIO

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

def filename_parsing(file):
    filename = file.split("/")[-1][:-4] # file name without extension
    sample_type = file.split("_")[1] # sample type: C - culture, P - plasmid etc.
    sample_name = file.split("_")[2] # sample name
    primer = file.split("_")[3] # primer
    return {'filename': filename, 'sample_type': sample_type, 'sample_name': sample_name, 'primer': primer, 'path': file}
    
def ab1_to_fastq(input_ab1, output_fq):
    record = SeqIO.read(input_ab1, "abi")
    SeqIO.write(record, output_fq, "fastq")

def run_bbduk(input_fq, output_fq, trimq, minlength):
    command = ['bbduk.sh',
           f'in={input_fq}',
           f'out={output_fq}',
           'qtrim=rl',
           f'trimq={trimq}',
           'qin=33',
           f'minlength={minlength}']
    p1 = subprocess.run(command, stderr=subprocess.DEVNULL)

#def fq_to_fa(input_fq, output_fa):

#def reverse_complement(input_fa, output_fa):

#def alignment(files):

#def find_consensus(aln):

def main():
    args = parse_arguments()
    #dir = args.input_directory
    #mode = args.files_processing_mode
    #if mode == "auto":

    dir = "../blast/data/101424"

    # Bulk processing
    files = []
    for item in os.listdir(dir):
        path=os.path.join(dir, item)
        if os.path.isfile(path):
            files.append(path)
    ab1_files=[file for file in files if ".ab1" in file] # full paths to ab1 files
    data = [filename_parsing(file) for file in ab1_files] # dictionary with files data

    # Pattern processing

    # Main cycle
    for file in data:
        # Generation of fastq files
        fq = file["path"][:-3]+"fq" # name of fastq file
        ab1_to_fastq(file["path"], fq) # create fq files from ab1 files

        # Trimming fastq files
        fq_trimmed = fq[:-3] + "_trimmed.fq" # name of trimmed fastq files
        run_bbduk(fq, fq_trimmed, trimq = 15, minlength = 50) # trimming fq files

    #path = "../blast/data/101424/Plate-2024-04-10_C_1_16SE1114-1096R_C02_03_2.ab1"

if __name__ == "__main__":
    main()