import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Reads trimming, build consensus and perform blastn.")
    parser.add_argument('--input_directory', type=str, default=".", help='Directory containing .ab1 files')
    parser.add_argument('--output_directory', type=str, default=".", help='Directory to save results')
    parser.add_argument('--files_processing_mode', type=str, default="auto", help='Files processing mode. Options: auto (default), manual')
    parser.add_argument('--filenames_patterns', type=list, default=list())
    parser.add_argument('--blastn_mode', type=str, default="auto", help='Balstn search mode. Options: auto (default), manual')
    parser.add_argument('--consensus_patterns', type=list, default=list())
    parser.add_argument('--nthreads', type=int, default=4, help='Number of threads for blastn search')
    parser.add_argument('--trimming_quality', type=int, default=15, help='Trimming quality')
    parser.add_argument('--minlength', type=int, default=50, help='Minimum length of consensus')
    parser.add_argument('--database', type=str, default=".", help='Blastn database')
    parser.add_argument('--qcovus', type=int, default=80, help='Query coverage threshold')
    parser.add_argument('--pident', type=int, default=95, help='Percent of identity threshold')
    return parser.parse_args()

def filenames_parsing(dir):
    args =  parse_arguments()
    dir = args.input_directory
    print(dir)
    return None

def main():

    tq = args.trimming_quality
    db = args.database
    print("dir is", dir)
    print("tq is", tq)
    print("db is", db)

if __name__ == "__main__":
    main()