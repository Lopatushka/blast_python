#!bin/bash
: '
Change substring in filenames.

Input: path to directory with files, pattern to be changed in filenames, substitute change.
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

