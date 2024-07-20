import os

def filename_parsing(file):
    sample_types = ["C", "CS", "F", "A"]
    #primer_names = []
    
    # Get filename without extension
    filename = file.split("/")[-1][:-4]

    # Split filename by _ and check its length
    filename_splitted = filename.split("_")
    if len(filename_splitted) != 7:
        raise ValueError(f"Error parsing filename {file}")
   
   # Get sample type and check it
    sample_type = filename_splitted[1]
    if sample_type not in sample_types:
        raise ValueError(f"Error parsing filename {file}: wrong sample_type: {sample_type}")
    
    # Get sample name
    sample_name = filename_splitted[2]

    # Get primer name
    primer = filename_splitted[3]
    primer_splitted = primer.split("-")

    # Check primer name
    if len(primer_splitted) != 2:
        raise ValueError(f"Error parsing filename {file}: wrong primer name: {primer}")
    #if primer_splitted[0] not in primer_names:
        #raise ValueError(f"Error parsing filename {file}: wrong primer name {primer}")
    if ("F" not in primer_splitted[1]) & ("f" not in primer_splitted[1]) & ("R" not in primer_splitted[1]) & ("r" not in primer_splitted[1]):
        raise ValueError(f"Error parsing filename {file}: wrong primer name: {primer}")

    # Get path to dir
    dir = "/".join(file.split("/")[:-1])
    # Check dir and parsing prosedure in general
    if not os.path.isdir(dir):
        raise ValueError(f"Error parsing filename {file}: directory {dir} doesn't exist")
    
    return {'filename': filename,
                'sample_type': sample_type,
                'sample_name': sample_name,
                'primer': primer,
                'path': file,
                'dir': dir}


def check_blast_database(database):
    dir = os.path.dirname(database)
    if not os.path.isdir(dir):
        raise IOError(f"Path {dir} to blastn database is not found")
    
    list_of_files = os.listdir(dir)

    # Check taxonomy_files are located in blast database dir
    taxonomy_files = ['taxdb.btd', 'taxdb.bti', 'taxonomy4blast.sqlite3']
    for taxonomy_file in taxonomy_files:
        if taxonomy_file not in list_of_files:
            raise IOError(f"{taxonomy_file} is not found in blastn database {database}")
    database_name = database.split("/")[-1]
    
    for taxonomy_file in taxonomy_files:
        list_of_files.remove(taxonomy_file)
    
    if sum([file[:-4] == database_name for file in list_of_files]) == 0:
        raise ValueError(f"Wrong value of blastn database name: {database_name}")


db = "../db_16S_ribosomal_RNA/16S_ribosomal_RNA"

check_blast_database(database=db)

    