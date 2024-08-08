#!bin/bash

FROM_DIR=$1

TARGET_DIR_NAME=$(basename "$FROM_DIR")
scp -r $FROM_DIR elopatuhina@10.100.20.211:/home/elopatuhina/sanger/${TARGET_DIR_NAME}
