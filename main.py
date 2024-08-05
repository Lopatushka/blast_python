from scripts import *

def main():
    print("Running program... Wait...")

    # Read command-line arguments
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

    # Check arguments provided - raising errors
    check_dir_exists(dir_path = dir)
    check_blast_database(database = database)
    check_percent_value(value_name = "--trimming_quality", value = trimming_quality)
    check_percent_value(value_name = "--qcovus", value = qcovus)
    check_percent_value(value_name = "--pident", value = pident)
    check_int_value(value_name = "--minlength", value = minlength, upper_bound = 500) # check if value < 0 or > 500
    check_int_value(value_name = "--nthreads", value = nthreads, upper_bound = 300) # check if value < 0 or > 300

    # Parsing starts
    if parsing_mode == "auto":
        ab1_files = list_of_files(dir, "ab1") # full paths to ab1 filesq
    elif parsing_mode == "manual":
        ab1_files = list_of_files_by_pattern(dir=dir, extension="ab1", patterns=parsing_patterns)
    
    data = [result for file in ab1_files if (result := filename_parsing(file))] # dictionary with files data

    # Initialize report dataframe
    report = pd.DataFrame(data)

    '''1st cycle - trimming, make revese complement'''
    for file in data:
        # Generation of fastq files
        fq = file["path"][:-3]+"fq" # path to fastq file
        ab1_to_fastq(file["path"], fq) # create fq file from ab1 file

        # Trimming fastq files
        fq_trimmed = fq[:-3] + "_trimmed.fq" # path to trimmed fastq file
        run_bbduk(fq, fq_trimmed, trimq = trimming_quality, minlength = minlength)

        # Convert trimmed fastq file to fasta file if exists AND not empty
        if os.path.isfile(fq_trimmed) and not_empty_file(fq_trimmed):
            # Add info to report
            report.loc[report["filename"] == file["filename"], "is_short"] = False
            fa_trimmed = fq_trimmed[:-2]+"fa" # path to trimmed fasta file
            fastq_to_fasta(fq_trimmed, fa_trimmed) # create trimmed fasta file

            # Make reverse complement if primer is reverse
            if ("R" in file['primer'].split("-")[1]) | ("r" in file['primer'].split("-")[1]):
            #if "R" in file['primer']:
                fa_trimmed_rc = fa_trimmed[:-3] + "_rc.fa" # path to fasta revese complement
                reverse_complement_fasta(fa_trimmed, fa_trimmed_rc) # make revesrse complement fasta file
                remove_file(fa_trimmed) # remove original fasta file
        else:
            # Add info to report
            report.loc[report["filename"] == file["filename"], "is_short"] = True
    
    remove_files_with_extension(dir = dir, extension="fq") # remove .fq files

    '''2nd cycle - make alignment, build consensus, blast'''
    if blastn_mode == "auto":
        fa_files = list_of_files(dir, "fa") # full paths to .fa files
    elif blastn_mode == "manual":
        fa_files = list_of_files_by_pattern(dir=dir, extension="fa", patterns=consensus_patterns)

    data = [result for file in fa_files if (result := filename_parsing(file))] # filename parsing of all .fa files
        
    sample_names = np.unique([file['sample_name'] for file in data if file]).tolist() #  store unique sample names in array
    
    # Store paths to .fa files with the identical sample names in array to build consensus if possible
    for sample_name in sample_names:
        pairs = [] 
        for file in data:
            if file['sample_name'] == sample_name:
                pairs.append(file['path'])
        length = len(pairs)

        if length == 1:
            # Add info to report
            report.loc[report["sample_name"] == sample_name, "is_consensus"] = False

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
            merge_name = f'{dir}/{sample_name}_merged.fa'
            merge_fasta_files(pairs, merge_name)

            # # Remove original files after merging to a single file
            # for item in pairs:
            #     remove_file(item)
            
            # Make alignment
            run_clustalw(merge_name)
            remove_files_with_extension(dir = dir, extension="dnd") # remove .dnd files

            # Make consensus
            aln_name = merge_name[:-2]+"aln" # path to aln file
            consensus_name = aln_name[:-4] + "_consensus.fa" # path to consensus file
            get_custom_consensus_from_aln(aln_name, consensus_name, threshold=0.7)

            # Check consensus quality
            if check_consensus_quality(consensus_name, threshold = consensus_quality): # good consensus
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
                
            else: # bad consensus
                # Add info to report
                report.loc[report["sample_name"] == sample_name, "is_consensus"] = False

                # Blastn search for files in pairs independently
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

    # Modify report dataframe
    report.drop(columns=["filename", "path", "dir"], inplace=True) # delete unnessesary cols
    report["is_consensus"] = report["is_consensus"].fillna(False)

    # Save report as .csv
    report.to_csv(dir + "/report.csv", sep='\t', index=False, header=True, encoding='utf-8', mode="a")

    print("Done!")


if __name__ == "__main__":
    main()