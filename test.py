import os
import subprocess
import io
import numpy as np
import pandas as pd
import warnings
from collections import Counter
from Bio import AlignIO
from Bio import SeqIO

def get_seqs_scores(fastq_files):
    seqs = []
    scores = []
    fastq_parser = SeqIO.parse(fastq_files, "fastq")
    for record in fastq_parser:
        seq = record.seq
        score = record.letter_annotations["phred_quality"]
        seqs.append(seq)
        scores.append(score)
    return(seqs, scores)

def get_custom_consensus_from_aln_score(aln_file, fq_files):
    try:
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
        seqs, scores = get_seqs_scores(fq_files)
        consensus = []
        
        for i in range(alignment_length):
            column = alignment[:, i]
            _count = column.count("-") # occurence of '-' in column
            column = column.replace("-", "") # delete "-" symbol from consensus     
            column += "-" *  _count # place '-' to the end of column
            counter = Counter(column) # create a speicial dict
            most_common_residue, count = counter.most_common(1)[0]
            if count / num_sequences >= threshold:
                consensus.append(most_common_residue)
            else:
                residue_set = frozenset([residue for residue in counter if residue != '-']) # delele '-' from column
                iupac_code = iupac_codes.get(residue_set, 'N')
                consensus.append(iupac_code)
        consensus_seq = ''.join(consensus) # make a string
        consensus_name = aln_file.split("/")[-1][:-4]
        # with open(consensus_fa, "w") as f: # save in fasta format
        #     f.write(f'>{consensus_name}_consensus' + "\n")
        #     f.write(consensus_seq)
    except ValueError as e:
        print(f"Error occured with file processing {aln_file}: {e}")
    except FileNotFoundError as e:
        print(f"Error occured with file processing {consensus_fa}: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")

fq1 = "./test/Plate-2024-04-10_C_1_16S155-F1_A02_01_2_trimmed.fq"
fq2 = "./test/Plate-2024-04-10_C_1_16SE1114.1096-R1_C02_03_2_trimmed.fq"
fq = "./test/fq_files.fq"



print(get_seqs_scores(fq))

