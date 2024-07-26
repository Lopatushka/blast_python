#!bin/bash

FROM_DIR=$1
TARGET_DIR_NAME=$2

scp -r $FROM_DIR elopatuhina@10.100.20.211:/home/elopatuhina/sanger_seq/${TARGET_DIR_NAME}
