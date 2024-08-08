#!/bin/bash

FOLDER_TO_SAVE=$1
LINK=$2

# Create the destination directory if it doesn't exist
mkdir -p $FOLDER_TO_SAVE

# Download the folder using gdown
gdown -O $FOLDER_TO_SAVE --folder $LINK

echo "Files sucessfully donwloaded and moved to folder $FOLDER_TO_SAVE"
