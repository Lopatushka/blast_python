import argparse
import os
import shutil
import subprocess
import sys
import io
import glob
import warnings
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
from collections import Counter

def parse_arguments():
    parser = argparse.ArgumentParser(description="Reads trimming, build consensus and perform blastn.")
    parser.add_argument('--directory', '-d', type=str, default=".", help='Directory containing .ab1 files')
    parser.add_argument('--parsing_mode', '-pm', type=str, default="auto", choices=["auto", "manual"], help='Filenames parsing mode. Options: auto (default), manual')
    parser.add_argument('--parsing_patterns', '-pp', nargs='*', type=str, default=list(), help='Unique patterns in filenames of .ab1 files for manual consensus building')
    parser.add_argument('--blastn_mode', '-bm', type=str, default="auto", choices=["auto", "manual"], help='Balstn search mode. Options: auto (default), manual')
    parser.add_argument('--consensus_patterns', '-cp', nargs='*', type=str, default=list(), help='Unique patterns in filenames of consensuses for manual blastn search')
    parser.add_argument('--consensus_quality', '-cq', type=int, default=15, help="Maximum percetage of 'N' in consensus")
    parser.add_argument('--nthreads', '-nt', type=int, default=4, help='Number of threads for blastn search')
    parser.add_argument('--trimming_quality', '-tq', type=int, default=15, help='Trimming quality')
    parser.add_argument('--minlength', '-ml', type=int, default=50, help='Minimum length of consensus')
    parser.add_argument('--database', '-db', type=str, default="/home/lopatushka/db/16S_ribosomal_RNA/16S_ribosomal_RNA", help='Blastn database')
    parser.add_argument('--qcovus', '-qc', type=int, default=80, help='Query coverage threshold')
    parser.add_argument('--pident', '-pi', type=int, default=95, help='Percent of identity threshold')
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
        
def merge_fasta_files(input_fasta_files, output_fasta):
    with open(output_fasta, "w") as output_handle:
        for input_file in input_fasta_files:
            records = SeqIO.parse(input_file, "fasta")
            SeqIO.write(records, output_handle, "fasta")

def run_clustalw(input_fasta):
    command = f'clustalw2 -INFILE={input_fasta}'
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
    with open(consensus_fa, "w") as f: # save in fasta format
        f.write(f'>{consensus_name}_consensus' + "\n")
        f.write(consensus_seq)

def check_consensus_quality(fasta, threshold):
    consensus = SeqIO.read(fasta, "fasta")
    consensus_length = len(consensus.seq)
    N_count = consensus.count("N")
    N_precentage = np.round(100*N_count/consensus_length, 3)
    if N_precentage >= threshold:
       return False # bad consensus
    else:
        return True # good consensus

def run_blastn(input_file, database, num_threads):
    command = f'blastn -query {input_file} -db {database} -num_threads {num_threads}\
        -outfmt "6 qseqid sacc staxid evalue pident mismatch gaps qcovus length sscinames"'
    p1 = subprocess.run(command, stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE, check=True, shell = True, text=True)
    return p1.stdout

def run_blastn_alignments(input_file, output_file, database, hits, num_threads):
    command = f'blastn -query {input_file} -db {database} -num_threads {num_threads} -taxids {hits} -outfmt 0 -out {output_file}'
    subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, shell = True, text=True)

def blastn_results_processing(data, qcovus_treshold, pident_treshold, consensus_name=None, database=None, dir="./"):
    # check axualiry files presence
    wd = os.getcwd() # working dir
    check_files_exist([wd + '/taxdb.btd', wd + '/taxdb.bti', wd + '/taxonomy4blast.sqlite3'])

    df = pd.read_csv(io.StringIO(data), index_col=False, header = None, sep = "\t",
                     names = ["qseqid", "sacc", "staxid", "evalue", "pident", "mismatch", "gaps", "qcovus", "length", "sscinames"])
    df["sscinames"] = df["sscinames"].str.split(" ").apply(lambda x: [str(x)] if isinstance(x, float) else x).apply(lambda x: x[:2]).apply(lambda x: " ".join(x))     
    # Filtration by query coverage - qcovus
    if sum(df["qcovus"] >= qcovus_treshold) < 5:
        df = df[:5] # show 5 first hits
    else:
        df = df[df["qcovus"] >= qcovus_treshold] # make subset

    # Insert sscinames_occurence column
    df.insert(9, "sscinames_occurence", 1)
    unique_values = df["sscinames"].value_counts()
    for inx in unique_values.index:
        df.loc[df["sscinames"] == inx, "sscinames_occurence"] = unique_values[inx]
    
    # Filtration by pident
    if sum(df["pident"] >= pident_treshold) < 5:
        df = df[:5] # show 5 first hits
    else:
        df = df[df["pident"] >= pident_treshold] # subset by pident
    
    df["pident"] = np.round(df["pident"], 3) # Round pident
    df = df.sort_values(by='pident', ascending=False) # Sort df by pident
    df = df[~(df["sscinames"].duplicated(keep='first'))] # Delete duplications of sscinames. Keep first occurences only
    
    # Get taxids of hits
    taxids = ','.join(map(str, list(df["staxid"])))

    # Save results
    output_file_tmp = dir + "/tmp.tsv"
    output_file = dir + "/blastn.tsv"
    df.to_csv(output_file_tmp, sep='\t', index=False, header=True, mode="w")
    with open(output_file_tmp, 'r') as file:
        content = file.read()
        add_line = f'file={consensus_name} db={database} qcovus_treshold={qcovus_treshold}, pident_treshold={pident_treshold}'
        content = add_line + "\n" + content + "\n"

    with open(output_file, 'a') as file: # save modified file
        file.write(content)
    
    if os.path.exists(output_file_tmp): # remove tmp file
        os.remove(output_file_tmp)
    
    return taxids

def is_empty(var, var_name):
    if not var:
        raise ValueError(f"Variable {var_name} is empty. Exiting the program.")

def check_1_file_exists(file_path):
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    
def check_files_exist(file_paths_list):
    for filepath in file_paths_list:
        check_1_file_exists(filepath)

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
        
def main():
    args = parse_arguments()
    dir = args.directory
    parsing_mode = args.parsing_mode
    parsing_patterns = args.parsing_patterns
    blastn_mode = args.blastn_mode
    consensus_patterns = args.consensus_patterns
    consensus_quality = args.consensus_quality
    nthreads = args.nthreads
    trimming_quality = args.trimming_quality
    minlength = args.minlength
    database = args.database
    qcovus = args.qcovus
    pident = args.pident
    
    #dir = "/home/lopatushka/blast/data/100424"
    #database = '/home/lopatushka/db/16S_ribosomal_RNA/16S_ribosomal_RNA'

    if parsing_mode == "auto":
        ab1_files = list_of_files(dir, "ab1") # full paths to ab1 files    
    elif parsing_mode == "manual":
        ab1_files = list_of_files_by_pattern(dir=dir, extension="ab1", patterns=parsing_patterns)
    
    is_empty(ab1_files, "ab1_files")

    data = [result for file in ab1_files if (result := filename_parsing(file))] # dictionary with files data
    is_empty(data, "data")

    # Add info to report dataframe
    report = pd.DataFrame(data)
         
    '''1st cycle - trimming, make revese complement'''
    for file in data:
        # Generation of fastq files
        fq = file["path"][:-3]+"fq" # path to fastq file
        ab1_to_fastq(file["path"], fq) # create fq file from ab1 file

        # Trimming fastq files
        fq_trimmed = fq[:-3] + "_trimmed.fq" # path to trimmed fastq file
        run_bbduk(fq, fq_trimmed, trimq = trimming_quality, minlength = minlength)

        # Convert trimmed fastq file to fasta file if exists
        if os.path.isfile(fq_trimmed):
            # Add info to report
            report.loc[report["filename"] == file["filename"], "is_short"] = False

            fa_trimmed = fq_trimmed[:-2]+"fa" # path to trimmed fasta file
            fastq_to_fasta(fq_trimmed, fa_trimmed) # create trimmed fasta file
            remove_file_by_pattern(file['dir'], pattern="*.fq") # remove .fq files

            # Make reverse complement if primer is reverse
            if "R" in file['primer']:
                fa_trimmed_rc = fa_trimmed[:-3] + "_rc.fa" # path to fasta revese complement
                reverse_complement_fasta(fa_trimmed, fa_trimmed_rc) # make revesrse complement fasta file
                remove_file(fa_trimmed) # remove original fasta file
        else:
            # Add info to report
            report.loc[report["filename"] == file["filename"], "is_short"] = True

    '''2nd cycle - make alignment, build consensus, blast'''
    if blastn_mode == "auto":
        fa_files = list_of_files(dir, "fa") # full paths to .fa files    
    elif blastn_mode == "manual":
        fa_files = list_of_files_by_pattern(dir=dir, extension="fa", patterns=consensus_patterns)
    
    is_empty(fa_files, "fa_files")
    
    data = [result for file in fa_files if (result := filename_parsing(file))] # filename parsing of all .fa files
    is_empty(data, "data")

    sample_names = np.unique([file['sample_name'] for file in data if file]).tolist() #  store unique sample names in array
    is_empty(sample_names, "sample_names")
    
    # Store paths to .fa files with the identical sample names in array to build consensus if possible
    for sample_name in sample_names:
        pairs = [] 
        for file in data:
            if file['sample_name'] == sample_name:
                pairs.append(file['path'])
        length = len(pairs)

        # Check how much files have the unique sample name: 0, 1 or more then one
        if length == 0: # in theory it is impossible situation
            warnings.warn(f"There is no files with the sample name {sample_name}")
            
        elif length == 1:
            # Blastn search for 1 file
            result = run_blastn(pairs[0], database, num_threads=nthreads)
            hits = blastn_results_processing(data=result, consensus_name=pairs[0],
                                            database=database, dir=dir,
                                            qcovus_treshold=qcovus, pident_treshold=pident)
            blast_aln = f'{dir}/blast_aln_{sample_name}.txt' # path to blastn_aln file
            run_blastn_alignments(input_file=pairs[0], output_file=blast_aln,
                                  database=database, hits=hits, num_threads=nthreads)
            
        else:
            # Merge files for alignment
            merge_name = f'{dir}/{sample_name}.fa'
            merge_fasta_files(pairs, merge_name)

            # Remove original files after merging to a single file
            for item in pairs:
                remove_file(item)
            
            # Make alignment
            run_clustalw(merge_name)
            remove_file_by_pattern(dir, '*.dnd') # delete .dnd files

            # Make consensus
            aln_name = merge_name[:-2]+"aln" # path to aln file
            consensus_name = aln_name[:-4] + "_consensus.fa" # path to consensus file
            get_custom_consensus_from_aln(aln_name, consensus_name, threshold=0.7)

            # Check consensus quality
            if check_consensus_quality(consensus_name, threshold = consensus_quality):
                # Add info to report
                report.loc[report["sample_name"] == sample_name, "is_consensus"] = True

                # Blastn search for consenus
                result = run_blastn(consensus_name, database, num_threads=nthreads)
                hits = blastn_results_processing(data=result, consensus_name=consensus_name,
                                                database=database, dir=dir,
                                                qcovus_treshold=qcovus, pident_treshold=pident)
                blast_aln = f'{dir}/blast_aln_{sample_name}.txt' # path to blastn_aln file
                run_blastn_alignments(input_file=consensus_name, output_file=blast_aln,
                                      database=database, hits=hits, num_threads=nthreads)
                
            else:
                # Add info to report
                report.loc[report["sample_name"] == sample_name, "is_consensus"] = False

                # Blastn search for files in pairs independently
                for file in pairs:
                    result = run_blastn(file, database, num_threads=nthreads)
                    hits = blastn_results_processing(data=result, consensus_name=consensus_name,
                                                    database=database, dir=dir,
                                                    qcovus_treshold=qcovus, pident_treshold=pident)
                    blast_aln = file[:-2] + 'txt' # path to blastn_aln file
                    run_blastn_alignments(input_file=file, output_file=blast_aln,
                                          database=database, hits=hits, num_threads=nthreads)
       
        pairs = []

    # Create new dirs and move files with special extensions to these dirs
    fa_files = list_of_files(dir, "fa") # full paths to fa files
    aln_files = list_of_files(dir, "aln") # full paths to aln files
    txt_files = list_of_files(dir, "txt") # full paths to txt files

    create_dir(dir + "/fasta")
    create_dir(dir + "/consensus_alns")
    create_dir(dir + "/blast_alns")

    move_files(dir + "/fasta", fa_files)
    move_files(dir + "/consensus_alns", aln_files)
    move_files(dir + "/blast_alns", txt_files)

    # Save report as .tsv
    report = report.drop(columns=["filename", "path", "dir"]) # delete unnessesary cols
    report.to_csv(dir + "/report.tsv", sep='\t', index=False, header=True, mode="a")


if __name__ == "__main__":
    main()