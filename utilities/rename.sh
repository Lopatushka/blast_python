#!bin/bash
: '
Description:
    This script is designed to rename files in a specified directory by changing a specific pattern in their filenames.
	The user provides a directory path, a pattern to search for in filenames, and the replacement string.
	The script scans the files in the given directory and modifies their filenames by substituting the specified pattern with the new string.

Usage:
    ./rename_files.sh <directory> <pattern> <substitute>

    <directory>   : Path to the directory containing files to be renamed.
    <pattern>     : Substring or pattern in the filenames that needs to be replaced.
    <substitute>  : The string that will replace the pattern in the filenames

Example:
    ./rename_files.sh /home/user/documents old-pattern new-pattern
'

DIR=$1
PATTERN=$2
RENAME=$3

FIND_RESULTS=$(find ${DIR} -maxdepth 1 -name "*${PATTERN}*")

while read FILE
do

	newfile=$(echo ${FILE} |sed "s/${PATTERN}/${RENAME}/")
	mv "${FILE}" "${newfile}" 

done <<< ${FIND_RESULTS}

