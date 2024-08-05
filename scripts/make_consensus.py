from Bio import SeqIO
from Bio import AlignIO
from collections import Counter
import subprocess
import numpy as np

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
    subprocess.run(command, stderr=subprocess.DEVNULL)

def fastq_to_fasta(input_fq, output_fa):
    with open(input_fq, "r") as fq, open(output_fa, "w") as fa:
        SeqIO.convert(fq, "fastq", fa, "fasta")

def reverse_complement_fasta(input_fasta, output_fasta):
    try:
        record = SeqIO.read(input_fasta, "fasta")
        record_rc = record.reverse_complement()
        record_rc.id = record.id
        record_rc.name = record.name
        with open(output_fasta, "w") as f:
            SeqIO.write(record_rc, f, "fasta")
    except ValueError as e:
        print(f"Error occured with file processing {input_fasta}: {e}")

def merge_fasta_files(input_fasta_files, output_fasta):
    try:
        with open(output_fasta, "w", encoding='utf-8') as output_handle:
            for input_file in input_fasta_files:
                records = SeqIO.parse(input_file, "fasta")
                SeqIO.write(records, output_handle, "fasta")
    except UnicodeEncodeError as e:
        print(f"Unicode encoding error during merging files: {input_fasta_files}. Error message: {e}")
    except Exception as e:
        print(f"An error occurred during merging files: {input_fasta_files}. Error message: {e}")

def run_clustalw(input_fasta):
    try:
        command = f'clustalw2 -INFILE={input_fasta}'
        subprocess.run(command, capture_output=True, check=True, shell = True)
    except subprocess.CalledProcessError as e:
        print(f"Error occured during clustalw run for file {input_fasta}. Error message: {e}")

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

def get_custom_consensus_from_aln(aln_file, consensus_fa, threshold=0.6):
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
        with open(consensus_fa, "w") as f: # save in fasta format
            f.write(f'>{consensus_name}_consensus' + "\n")
            f.write(consensus_seq)
    except ValueError as e:
        print(f"Error occured with file processing {aln_file}: {e}")
    except FileNotFoundError as e:
        print(f"Error occured with file processing {consensus_fa}: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")

def check_consensus_quality(fasta, threshold):
    try:
        consensus = SeqIO.read(fasta, "fasta")
        consensus_length = len(consensus.seq)
        N_count = consensus.count('R') + consensus.count('Y') + consensus.count('S')\
                            + consensus.count('W') + consensus.count('K') + consensus.count('M')\
                            + consensus.count('V') + consensus.count('H') + consensus.count('D')\
                            + consensus.count('B') + consensus.count("N")
        N_precentage = np.round(100*N_count/consensus_length, 3)
        if N_precentage >= threshold:
            return False # bad consensus
        else:
            return True # good consensus
    except FileNotFoundError as e:
        print(f"Error occured with file processing {fasta}: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")

