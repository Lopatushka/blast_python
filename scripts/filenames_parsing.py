import argparse
import os
import subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
from collections import Counter
import io

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
    with open(output_fasta, "w") as f:
            SeqIO.write(record_rc, f, "fasta")

def remove_file_by_pattern(dir, pattern):
    command = f'rm {dir}/*{pattern}*'
    subprocess.run(command, capture_output=True, check=True, shell=True)

def remove_file(path_to_file):
    command = f'rm -f {path_to_file}'
    subprocess.run(command, capture_output=True, check=True, shell=True)
        
def merge_fasta_files(input_fasta_files, output_fasta):
    with open(output_fasta, "w") as output_handle:
        for input_file in input_fasta_files:
            records = SeqIO.parse(input_file, "fasta")
            SeqIO.write(records, output_handle, "fasta")

def run_clustalw(input_fasta):
    command = f'clustalw2 {input_fasta}'
    subprocess.run(command, capture_output=True, check=True, shell = True)

def get_custom_consensus_from_aln(aln_file, consensus_fa, threshold=0.7):
    # IUPAC nucleotide codes
    iupac_codes = {
    frozenset(['A']): 'A',
    frozenset(['C']): 'C',
    frozenset(['G']): 'G',
    frozenset(['T']): 'T',
    frozenset(['A', 'G']): 'R',
    frozenset(['C', 'T']): 'Y',
    frozenset(['G', 'C']): 'S',
    frozenset(['A', 'T']): 'W',
    frozenset(['G', 'T']): 'K',
    frozenset(['A', 'C']): 'M',
    frozenset(['A', 'C', 'G']): 'V',
    frozenset(['A', 'C', 'T']): 'H',
    frozenset(['A', 'G', 'T']): 'D',
    frozenset(['C', 'G', 'T']): 'B',
    frozenset(['A', 'C', 'G', 'T']): 'N',
}
    alignment = AlignIO.read(aln_file, "clustal")
    num_sequences = len(alignment) # number of sequences
    alignment_length = alignment.get_alignment_length() # length of each sequences
    consensus = []
    for i in range(alignment_length):
        column = alignment[:, i]
        _count = column.count("-") # occurence of '-' in column
        column = column.replace("-", "") # delete "-" symbol from consensus     
        column += "-" *  _count # place '-' to the end of column
        counter = Counter(column)
        most_common_residue, count = counter.most_common(1)[0]
        if count / num_sequences >= threshold:
            consensus.append(most_common_residue)
        else:
            residue_set = frozenset([residue for residue in counter if residue != '-']) # delele '-' from column
            iupac_code = iupac_codes.get(residue_set, 'N')
            consensus.append(iupac_code)
    consensus_seq = ''.join(consensus) # make a string
    consensus_name = aln_file.split("/")[-1][:-4]
    with open(consensus_fa, "w") as f:
        f.write(f'>{consensus_name}_consensus' + "\n")
        f.write(consensus_seq)

def check_consensus_quality(fasta, threshold = 15):
    consensus = SeqIO.read(fasta, "fasta")
    consensus_length = len(consensus.seq)
    N_count = consensus.count("N")
    N_precentage = np.round(100*N_count/consensus_length, 3)
    if N_precentage >= threshold:
       return False # bad consensus
    else:
        return True # good consensus

def run_blastn(input_file, database, num_threads=4):
    command = f'blastn -query {input_file} -db {database} -num_threads {num_threads}\
        -outfmt "6 qseqid sacc staxid evalue pident mismatch gaps qcovus length sscinames"'
    p1 = subprocess.run(command, stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE, check=True, shell = True, text=True)
    return p1.stdout 

#def blastn_results_processing(data):

def main():
    args = parse_arguments()
    #dir = args.input_directory
    #mode = args.files_processing_mode
    #if mode == "auto":

    dir = "../blast/data/101424"
    database = '../db_16S_ribosomal_RNA/16S_ribosomal_RNA'

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
            # Merge files for alignment
            subset_pairs = [filename_parsing(file) for file in pairs]
            subset_sample_name = np.unique([file['sample_name'] for file in subset_pairs])[0]
            subset_sample_dir = np.unique([file['dir'] for file in subset_pairs])[0]
            merge_name = f'{subset_sample_dir}/{subset_sample_name}.fa'
            merge_fasta_files(pairs, merge_name)
            
            # Make alignment
            run_clustalw(merge_name)
            remove_file_by_pattern(subset_sample_dir, 'dnd') # delete .dnd files

            # Make consensus
            aln_name = merge_name[:-2]+"aln" # path to aln file
            consensus_name = aln_name[:-4] + "_consensus.fa" # path to consensus file
            get_custom_consensus_from_aln(aln_name, consensus_name, threshold=0.7)

            # Check consensus quality
            if check_consensus_quality(consensus_name, threshold = 15):
                # Blastn search for consenus
                result = run_blastn(consensus_name, database, num_threads=4)
                
            else:
                # Blastn search for files in pairs independently
                pass

            #print(pairs, length)

        
        pairs = []



if __name__ == "__main__":
    main()