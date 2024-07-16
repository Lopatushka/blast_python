```markdown
# Automatic Sanger sequencing results analysis

## Description
This package makes automated bulk Sanger sequencing data analysis.

## Prerequisites
python >=3.8.19 version or higher
pandas 2.0.3
numpy 1.24.4
biopython 1.83
bbduk.sh version 39.08
CLUSTAL 2.1
BLAST 2.15.0+

## Instalation

## Algorithm
### Trimming Sanger sequencing reads
Trimming is performed using BBMap software from both right and left ends by Phred algotithm with following arguments:

--trimming_quality specifies trimq argument  (default is 15)
--minlength specifies minlength argument (default is 50)

For details see BBMap manual.

### Consensus finding
Consensus is built using custom defined funtion.
If the number of 'N' in consensus file exceeds consensus_quality threshold (--consensus_quality) (default = 15%), consensus is qualified is bad.
In this case, blastn search will be performed for F and R reads independently.

### Blastn search
Blastn search is performed against locally installed blast database. You need to specify the path to this database (--database)
Resalts are stored in ./blastn.tsv table.
Results are sorted by pident and filtered subsequently by query coverage (--qcovus, default = 80%) and percent of identity (--pident, default = 95%).
If there is no hits after filtration, the program shows first five hits.

For details see blastn manual.

## Usage
To use the package launch the program from the commande line.
Example of code for linux-command line:

```bash
dir = './tests'
database = '/home/lopatushka/db/16S_ribosomal_RNA/16S_ribosomal_RNA'
python -m blast_python.main -d $dir -db $database

## Nucleotide Blast database
You can use a ready-made NCBI blast database or build your own database.
To check all blast databases allowed for downloading use the following command:
```bash
update_blastdb.pl --showal

To download database use this:
```bash
mkdir -p blast_database cd blast_database update_blastdb.pl name_of_database --decompress

If you need to use taxonomy-related information:
```bash
update_blastdb.pl taxdb --decomress

To build your own blast database with taxonomy information you need to follow these steps:
1. Download sequences of interest in fasta format to a single .fa file.
2. Make a special file which contains information about accession numbers and taxonomy ids.
3. Create blast database
4. Download taxdb.bti, taxdb.btd, and taxonomy4blast.sqlite3 files and put them to the same directory where the blast database is located.
5. Check your blast database

For details see BLAST+ manual.

### 1. Fasta sequence downloading
You can use both web version of NCBI Nucleotide database or use E-Utilities.

Example:
```bash
mkdir -p my_database cd ./my_database
esearch -db nuccore -query "16S[All Fields] AND rRNA[All Fields] NOT Uncultured[All Fields] \
                            AND(feces[All Fields] OR stool[All Fields] OR gut[All Fields] OR fecal[All Fields]) \
                            AND ("1000"[SLEN]:"2000"[SLEN])" \
                            | efetch -format fasta > files.fa

For details see E-Utilities manual, esearch -help.

### 2. Taxonomy mask file
To add information about taxonomy to blast database you need to create a special file with two columns:
first column is the sequence ID, the second column should be the NCBI taxonomy ID.
You can make this file manually or download a special file called nucl_gb.accession2taxid.gz from here:
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/

First, download the file as .gz archive.
Next, decompress the file and check it.
Finally, exctract two columns you needed from this file, delete the first string and write results in a new file.

Example of code:
```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz gunzip
cat nucl_gb.accession2taxid | head -5
cut -f 2,3 nucl_gb.accession2taxid | awk 'NR>1' > nucl_gb.accession2taxid_use

Example of mask file with two rows:
LN998086.1 10641
AB117566.1 9913

For details see BLAST+ manual and makeblastdb manual in particularly:
https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.makeblastdb_application_opt/

### 3. Build blast indexes
To create blast database use the command makeblastdb.
Example:
```bash
makeblastdb -in files.fa -parse_seqids -taxid_map mask_file -dbtype nucl -out my_database

For details see: makeblastdb -help, makeblastdb manual:
https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.makeblastdb_application_opt/

### 4. taxdb.bti, taxdb.btd, and taxonomy4blast.sqlite3 files
Download taxdb.bti, taxdb.btd, and taxonomy4blast.sqlite3 files and put them to the directory of your database:
```bash
update_blastdb.pl taxdb --decomress

For details see BLAST+ manual.

### 5. Check blast database
Use the following commands for blast database checking:
```bash
blastdbcmd -db my_databse -info blastdbcmd -db my_databse -tax_info

For details see: BLAST+ manual, blastdbcmd -help






