import os
import subprocess
import io
import numpy as np
import pandas as pd
import warnings
from collections import Counter
from Bio import AlignIO
from Bio import SeqIO

def check_consensus_quality(fasta, threshold):
    try:
        consensus = SeqIO.read(fasta, "fasta")
        consensus_length = len(consensus.seq)
        N_count = consensus.count('R') + consensus.count('Y') + consensus.count('S')\
              + consensus.count('W') + consensus.count('K') + consensus.count('M')\
                + consensus.count('V') + consensus.count('H') + consensus.count('D') + consensus.count('B') + consensus.count("N")
        N_precentage = np.round(100*N_count/consensus_length, 3)
        print(N_precentage)
        if N_precentage >= threshold:
            return False # bad consensus
        else:
            return True # good consensus
    except FileNotFoundError as e:
        print(f"Error occured with file processing {fasta}: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")

path = "../sanger_seq/310724/fasta/13A_consensus.fa"
threshold = 10

s = '../sanger_seq/310724/Plate-2024-07-31_CS_7A_strep.sp-R2_D02_04_trimmed_rc.fa'
print(os.path.basename(s))

