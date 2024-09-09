#!bin/bash
: '
Description:
    This script is designed to copy a directory from a source path (within the Windows Subsystem for Linux (WSL) environment)
    to a target directory on the local Linux filesystem.
    The directory to be copied is specified by the user as a command-line argument,
    and the script creates the target directory if it does not exist before performing the copy operation.
'

DIR_NAME=$1
#TARGET_DIR=$2

FROM_DIR=/mnt/c/users/ELopatuhina/Downloads/${DIR_NAME}
TARGET_DIR=/home/elopatuhina/sanger/${DIR_NAME}

mkdir -p "${TARGET_DIR}"

cp -r $FROM_DIR $TARGET_DIR

