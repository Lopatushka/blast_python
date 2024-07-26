import os
import subprocess
import io
import numpy as np
import pandas as pd
import warnings
from .errors_catching import check_files_exist

def run_blastn(input_file, database, num_threads):
    try:
        command = f'blastn -query {input_file} -db {database} -num_threads {num_threads}\
            -outfmt "6 qseqid sacc staxid evalue pident mismatch gaps qcovus length sscinames"'
        p1 = subprocess.run(command, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, check=True, shell = True, text=True)
        return p1.stdout
    except subprocess.CalledProcessError as e:
        print(f"There is no match for query {input_file} in blastn database {database}")
        print(f"Error message is: {e}")

def run_blastn_alignments(input_file, output_file, database, hits, num_threads):
    try:
        command = f'blastn -query {input_file} -db {database} -num_threads {num_threads} -taxids {hits} -outfmt 0 -out {output_file}'
        subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, shell = True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occured during blastn alignments run for query {input_file} and blastn database {database}")
        print(f"Error message is: {e}")

def blastn_results_processing(data, qcovus_treshold, pident_treshold, consensus_name=None, database=None, dir="./"):
    try:
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
        output_file_tmp = dir + "/tmp.csv"
        output_file = dir + "/blastn.csv"
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
    except Exception as e:
        print(f"Unexpected erroe occured: {e}")