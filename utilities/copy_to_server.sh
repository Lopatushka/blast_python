#!bin/bash
: '
    This script is designed to securely copy (using SCP) a directory from the local Windows Subsystem for Linux (WSL) environment to a remote server.
    The directory to be copied is specified by the user as a command-line argument,
    and the script transfers it to a specified directory on the remote server using the `scp` command.

'

DIR_NAME=$1
FROM_DIR=/mnt/c/users/ELopatuhina/Downloads/${DIR_NAME}

echo "Copy directory ${DIR_NAME} to server"

scp -r $FROM_DIR elopatuhina@10.100.20.211:/home/elopatuhina/sanger/${DIR_NAME}
