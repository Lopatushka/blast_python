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

def list_of_files(dir, extension):
        files = []
        for item in os.listdir(dir):
            path=os.path.join(dir, item)
            if os.path.isfile(path):
                files.append(path)
        subset_of_files=[file for file in files if f'.{extension}' in file] # full paths to .ext files
        return subset_of_files


def filename_parsing(file):
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

def fastq_to_fasta(input_fq, output_fa):
    with open(input_fq, "r") as fq, open(output_fa, "w") as fa:
        SeqIO.convert(fq, "fastq", fa, "fasta")

def reverse_complement_fasta(input_fasta, output_fasta):
    record = SeqIO.read(input_fasta, "fasta")
    record_rc = record.reverse_complement()
    record_rc.id = record.id
    record_rc.name = record.name
    with open(output_fasta, "a") as f:
            SeqIO.write(record_rc, f, "fasta")

def remove_file_by_pattern(dir, pattern):
    command = f'rm {dir}/*{pattern}*'
    subprocess.run(command, capture_output=True, check=True, shell=True)

def remove_file(path_to_file):
    command = f'rm -f {path_to_file}'
    subprocess.run(command, capture_output=True, check=True, shell=True)
        
#def run_clustalw(files):

#def find_consensus(aln):

def main():
    args = parse_arguments()
    #dir = args.input_directory
    #mode = args.files_processing_mode
    #if mode == "auto":

    dir = "../blast/data/101424"

    # Bulk processing
    ab1_files = list_of_files(dir, "ab1") # full paths to ab1 files
    data = [filename_parsing(file) for file in ab1_files] # dictionary with files data

    # Processing by pattern
    # ab1_files =
    # data = 

    # 1st cycle - trimming, make revese complement
    for file in data:
        #print(file)
        # Generation of fastq files
        fq = file["path"][:-3]+"fq" # path to fastq file
        ab1_to_fastq(file["path"], fq) # create fq file from ab1 file

        # Trimming fastq files
        fq_trimmed = fq[:-3] + "_trimmed.fq" # path to trimmed fastq file
        run_bbduk(fq, fq_trimmed, trimq = 15, minlength = 50)

        # Convert trimmed fastq file to fasta file
        fa_trimmed = fq_trimmed[:-2]+"fa" # path to trimmed fasta file
        fastq_to_fasta(fq_trimmed, fa_trimmed) # create trimmed fasta file
        remove_file_by_pattern(file['dir'], pattern=".fq") # remove .fq files

        # Make reverse complement if primer is reverse
        if "R" in file['primer']:
            fa_trimmed_rc = fa_trimmed[:-3] + "_rc.fa" # path to fasta revese complement
            reverse_complement_fasta(fa_trimmed, fa_trimmed_rc) # make revesrse complement fasta file
            remove_file(fa_trimmed) # remove original fasta file

    # 2nd cycle - make alignment, build consensus
    fa_files = list_of_files(dir, "fa") # full paths to .fa files
    data = [filename_parsing(file) for file in fa_files]
    sample_names = np.unique([file['sample_name'] for file in data])
    for sample_name in sample_names:
        pairs = []
        for file in data:
            if file['sample_name'] == sample_name:
                pairs.append(file['path'])
        length = len(pairs)
    
        if length == 1:
            # blast
            pass
        else:
            # make aln
            # make consensus
            # blast
            pass

        print(pairs, length)
        pairs = []


if __name__ == "__main__":
    main()