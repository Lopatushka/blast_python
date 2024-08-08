#!bin/bash

DIR_NAME=$1

#FROM_DIR=/mnt/c/users/ELopatuhina/Desktop/Work/Data/Raw_data/seq/${DIR_NAME}
FROM_DIR=/mnt/c/users/ELopatuhina/Downloads/${DIR_NAME}
TARGET_DIR=/home/lopatushka/sanger_seq

mkdir -p "${TARGET_DIR}"

cp -r $FROM_DIR $TARGET_DIR

