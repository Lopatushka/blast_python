from scripts import *

def main():
    # Initialize arguments
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

    print("Running program... Wait...")

    #print(args)
    
    # Parsing starts
    if parsing_mode == "auto":
        ab1_files = list_of_files(dir, "ab1") # full paths to ab1 filesq
    elif parsing_mode == "manual":
        ab1_files = list_of_files_by_pattern(dir=dir, extension="ab1", patterns=parsing_patterns)

    # Check if there is ab1 files founded
    if not ab1_files:
        raise FileNotFoundError(f"No .ab1 files is founded by your query.")
    
    # Dictionary with ab1 files metadata
    data = [result for file in ab1_files if (result := filename_parsing(file))]
    
    # Initialize report dataframe, add two additional columns
    report = pd.DataFrame(data)
    report["is_long"] = False
    report["is_consensus"] = False
    
    '''1st cycle - trimming, make revese complement'''
    for file in data:
        # Generation of fastq files
        fq = file["path"][:-3]+"fq" # path to fastq file
        ab1_to_fastq(file["path"], fq) # create fq file from ab1 file

        # Trimming fastq files
        fq_trimmed = fq[:-3] + "_trimmed.fq" # path to trimmed fastq file
        run_bbduk(fq, fq_trimmed, trimq = trimming_quality, minlength = minlength)

        # Convert trimmed FASTQ file to FASTA file if it exists AND not empty
        if is_file_exists(fq_trimmed):
            # Add info to report
            report.loc[report["filename"] == file["filename"], "is_long"] = True
            fa_trimmed = fq_trimmed[:-2]+"fa" # path to trimmed fasta file
            fastq_to_fasta(fq_trimmed, fa_trimmed) # create trimmed fasta file

            # Make reverse complement if primer is reverse
            if ("R" in file['primer_orientation']) | ("r" in file['primer_orientation']):
                fa_trimmed_rc = fa_trimmed[:-3] + "_rc.fa" # path to fasta revese complement
                reverse_complement_fasta(fa_trimmed, fa_trimmed_rc) # make revesrse complement fasta file
                remove_file(fa_trimmed) # remove original fasta file
    
    # Remove all FASTQ files
    remove_files_with_extension(dir = dir, extension="fq")

    '''2nd cycle - make alignment, build consensus, blast'''
    if blastn_mode == "auto":
        fa_files = list_of_files(dir, "fa") # full paths to .fa files
        if not fa_files:
            warnings.warn("There is no FASTA files for BLASTN search. Check --minlength value provided.")
            return None
    elif blastn_mode == "manual":
        fa_files = list_of_files_by_pattern(dir=dir, extension="fa", patterns=consensus_patterns)
        if not fa_files:
            warnings.warn("There is no FASTA files for BLASTN search. Check --consensus_patterns arguments and --minlength value provided.")
            return None
        
    # Dictionary with FASTA files metadata   
    data = [result for file in fa_files if (result := filename_parsing(file))]
    # Store unique sample_names_primers in array
    sample_names = np.unique([file['sample_name_primer'] for file in data if file]).tolist()
    
    # Store paths to FASTA files with the identical sample names in array to build consensus if possible
    for sample_name in sample_names:
        pairs = [] 
        for file in data:
            if file['sample_name_primer'] == sample_name:
                pairs.append(file['path'])
        length = len(pairs)

        # ONLY ONE FASTA file with this particular 'sample_name_primer'
        if length == 1:
            # Add info to report
            #report.loc[report["sample_name_primer"] == sample_name, "is_consensus"] = False

            # BLASTN search
            result = run_blastn(pairs[0], database, num_threads=nthreads)
            hits = blastn_results_processing(data=result, consensus_name=pairs[0],
                                            database=database, dir=dir,
                                            qcovus_treshold=qcovus, pident_treshold=pident)
            blast_aln = f'{dir}/blast_aln_{sample_name}.txt' # path to blastn_aln file
            run_blastn_alignments(input_file=pairs[0], output_file=blast_aln,
                                  database=database, hits=hits, num_threads=nthreads)
            
        # MORE THEN ONE FASTA file with this particular 'sample_name_primer'
        else:
            # Merge files for alignment
            merge_name = f'{dir}/{sample_name}.fa'
            merge_fasta_files(pairs, merge_name)

            # Make alignement
            run_clustalw(merge_name)
            remove_files_with_extension(dir = dir, extension="dnd") # remove .dnd files

            # Make consensus
            aln_name = merge_name[:-2]+"aln" # path to aln file
            consensus_name = aln_name[:-4] + "_consensus.fa" # path to consensus file
            get_custom_consensus_from_aln(aln_name, consensus_name, threshold=0.6)

            # Check consensus quality
            # If GOOD consensus
            if check_consensus_quality(consensus_name, threshold = consensus_quality): 
                # Add info to report
                report.loc[report["sample_name_primer"] == sample_name, "is_consensus"] = True

                # BLASTN search for consenus
                result = run_blastn(consensus_name, database, num_threads=nthreads)
                hits = blastn_results_processing(data=result, consensus_name=consensus_name,
                                                database=database, dir=dir,
                                                qcovus_treshold=qcovus, pident_treshold=pident)
                blast_aln = f'{dir}/blast_aln_{sample_name}.txt' # path to blastn_aln file
                run_blastn_alignments(input_file=consensus_name, output_file=blast_aln,
                                      database=database, hits=hits, num_threads=nthreads)
            # If BAD consensus
            else:
            #     # Add info to report
            #     report.loc[report["sample_name_primer"] == sample_name, "is_consensus"] = False

                # BLASTN search for files in pairs independently
                for file in pairs:
                    result = run_blastn(file, database, num_threads=nthreads)
                    hits = blastn_results_processing(data=result, consensus_name=os.path.basename(file),
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

    # Delete unnessesary cols in report
    report.drop(columns=["filename", "path", "dir", "primer_orientation"], inplace=True)
    # report["is_consensus"] = report["is_consensus"].fillna(False)

    # Save report as .csv
    report.to_csv(dir + "/report.csv", sep='\t', index=False, header=True, encoding='utf-8', mode="a")

    print("Done!")


if __name__ == "__main__":
    main()