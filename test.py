import os
import subprocess
import io
import numpy as np
import pandas as pd
import warnings
from collections import Counter
from Bio import AlignIO
from Bio import SeqIO

def run_blastn(input_file, database, num_threads):
    try:
        command = f'blastn -query {input_file} -db {database} -num_threads {num_threads}\
            -outfmt "6 qseqid sacc staxid evalue pident mismatch gaps qcovus length sscinames"'
        p1 = subprocess.run(command, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, check=True, shell = True, text=True)
        if not p1.stdout:
            print(f"There is no match for query {input_file} in blastn database {database}")
        return p1.stdout
    except subprocess.CalledProcessError as e:
        print(f"There is no match for query {input_file} in blastn database {database}")

print(run_blastn("./test/test.fa" , "../db/16S_ribosomal_RNA/16S_ribosomal_RNA", 4))