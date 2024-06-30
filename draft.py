#!/usr/bin/env python3
import os
import subprocess
import sys
import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align.AlignInfo import SummaryInfo
from Bio import AlignIO
import pandas as pd
import numpy as np
from Bio.Blast.Applications import NcbiblastnCommandline
import io 

def reverse_complement(args, path_to_save):
    num_args = len(args)
    if num_args == 2:
        record1=SeqIO.read(args[0], "fastq")
        record2=SeqIO.read(args[1], "fastq")
        record_rc=record2.reverse_complement()
        record_rc.id = record2.id
        record_rc.name = record2.name
        record_rc.description=record2.description
        records=[record1, record_rc]
        with open(path_to_save, "a") as handle:
            SeqIO.write(records, handle, "fasta")

    elif num_args == 1:
        record=SeqIO.read(args[0], "fastq")
        with open(path_to_save, "a") as handle:
            SeqIO.write(record, handle, "fasta")

def find_consensus(aln, threshold, path_to_save):
    align = AlignIO.read(aln, "clustal")
    consensus = SummaryInfo(align).dumb_consensus(0.7, "N") # find consensus
    consensus = str(consensus)
    consensus_length = len(consensus)
    N_count = consensus.count("N")
    N_precentage = np.round(100*N_count/consensus_length, 3)

    if N_precentage >= threshold:
       return False
    else:
        consensus_name = ">" + aln.split("/")[-1].rsplit(".")[0] + "_Consensus"
        lines = [consensus_name, consensus]
        consensus_fasta = "\n".join(lines)
        with open(path_to_save, "w") as handle:
            handle.write(consensus_fasta)
        return True

#def trimming_consensus():


def run_blastn(dir, db, qcovus_treshold, pident_treshold):
    # Get all files from the dir
    dir_to_files = dir + "/blast_files"
    files = []
    for item in os.listdir(dir_to_files):
        path=os.path.join(dir_to_files, item)
        if os.path.isfile(path):
            files.append(path)

    fasta_files=[file for file in files if ".fa" in file] 

    for file in fasta_files:
        blastn_cline = NcbiblastnCommandline(query = file, db = db, outfmt = "6 qseqid sacc staxid pident mismatch gaps qcovus length sscinames")
        stdout_tbl, stderr_tbl = blastn_cline()
        df = pd.read_csv(io.StringIO(stdout_tbl), index_col=False, header = None, sep = "\t", names = ["qseqid", "sacc", "staxid", "pident", "mismatch", "gaps", "qcovus", "length", "sscinames"])

        df["sscinames"] = df["sscinames"].str.split(" ").apply(lambda x: [str(x)] if isinstance(x, float) else x).apply(lambda x: x[:2]).apply(lambda x: " ".join(x))     

        # Filtration by query coverage - qcovus
        if sum(df["qcovus"] >= qcovus_treshold) < 5:
            df = df[:5] # show 5 first hits
        else:
            df = df[df["qcovus"] >= qcovus_treshold] 

        # Insert sscinames_occurence column
        df.insert(8, "sscinames_occurence", 1)
        unique_values = df["sscinames"].value_counts()
        for inx in unique_values.index:
            df.loc[df["sscinames"] == inx, "sscinames_occurence"] = unique_values[inx]

        # Filtration by pident
        if sum(df["pident"] >= pident_treshold) < 5:
            df = df[:5] # show 5 first hits
        else:
            df = df[df["pident"] >= pident_treshold] 
    
        # Round pident
        df["pident"] = np.round(df["pident"], 3)

        # Delete duplications of sscinames. Keep first occurences only
        df = df[~(df["sscinames"].duplicated(keep='first'))]

        taxids = ','.join(map(str, list(df["staxid"])))
        #blastn_cline = NcbiblastnCommandline(query = file, db = db, outfmt = "0", taxids = taxids)
        #stdout_aln, stderr_aln = blastn_cline()
             
        df.to_csv(dir + "/blast_tbl.txt", sep='\t', index=False, header=False, mode='a')# make tabular separated string


def main():
    dir = sys.argv[1]

    # Get all files from the dir
    files = []
    for item in os.listdir(dir):
        path=os.path.join(dir, item)
        if os.path.isfile(path):
            files.append(path)

    ab1_files=[file for file in files if ".ab1" in file]
    sample_names=np.unique([file.split("_")[1] for file in ab1_files])

    # Make pairs of files
    array=[]
    for name in sample_names:
        args=[]
        for file in ab1_files:
            if name == file.split("_")[1]:
                args.append(file)
                continue
        array.append(args)    

    # File processing - trimming
    for args in array:
        length = len(args)
        print("Number of files to process:", length)
        print("Files to process:", *args)

        filenames =  [file.split("/")[-1] for file in args]
        sample_name = "".join(np.unique([filename.split("_")[1] for filename in filenames]))
        paths_to_fq = [dir + "/blast_files/" + name[:-4] + ".fq" for name in filenames]
        paths_to_fq_trimmed = [dir + "/blast_files/" + name[:-4] + "_trimmed.fq" for name in filenames]
        path_to_fa= dir + "/" + sample_name + ".fa"

        print("Creating blast_files directory")
        os.system("mkdir -p" + " " +  dir + "/blast_files")
       
        for (file, fq, fq_trimmed) in zip(args, paths_to_fq, paths_to_fq_trimmed):
            os.system("seqret -sformat abi -osformat fastq -auto -stdout -sequence" + " " + file + ">" + fq)
            os.system("bbduk.sh -Xmx2g in=" + fq + " " + "out=" + fq_trimmed + " qtrim=rl trimq=15 qin=33 minlength=50 > /dev/null 2>&1")
        
        fq_trimmed_bool = [os.path.exists(fq_trimmed) for fq_trimmed in paths_to_fq_trimmed]
        fq_trimmed_existed = [x for x, y in zip(paths_to_fq_trimmed, fq_trimmed_bool) if y]
        fq_trimmed_existed_names = [name.split("/")[-1][:-11] for name in fq_trimmed_existed]
        n_fq_trimmed = len(fq_trimmed_existed)

        # Convert fastq files to fasta
        #for file in fq_trimmed_existed:
         #   SeqIO.convert(file, "fastq", file.replace("_trimmed.fq", ".fa"), "fasta")
        
        if n_fq_trimmed == 2: # make parsing and > 1 condition
            print("Both reads are high quality. Make alignement and build consensus.")
            reverse_complement(fq_trimmed_existed, path_to_fa) # convert to fasta and make reverse complement
            os.system("clustalw" + " " + path_to_fa + " -QUIET > /dev/null 2>&1")
            path_to_aln = dir + "/" + sample_name + ".aln"
            path_to_consensus = dir + "/blast_files/" + sample_name + ".fa"
            if find_consensus(path_to_aln, 10, path_to_consensus): # build consensus
                os.system("rm " + dir + "/*dnd")
                os.system("rm " + dir + "/*fa")
                
                print("Creating aligns_for_consensus directory")
                os.system("mkdir -p " + dir + "/aligns_for_consensus")
                os.system("mv " + dir + "/*aln " + dir + "/aligns_for_consensus")

        elif n_fq_trimmed == 0:
            print("All reads are low quality!")

        elif n_fq_trimmed == 1:
            print("No consensus building. Blast search for read: ", *fq_trimmed_existed_names)
            reverse_complement(fq_trimmed_existed, path_to_fa)

    os.system("rm " + dir + "/blast_files/*fq")

    db = "../db_16S_ribosomal_RNA/16S_ribosomal_RNA"
    run_blastn(dir, db, 80, 95)

if __name__ == "__main__":
    main()
