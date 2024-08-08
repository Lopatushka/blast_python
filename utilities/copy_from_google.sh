#!/bin/bash

FOLDER_ID=$1
FOLDER_TO_SAVE=$2

# Create the destination directory if it doesn't exist
mkdir -p $FOLDER_TO_SAVE

# Download the folder using gdown
gdown -O $FOLDER_TO_SAVE --folder https://drive.google.com/drive/folders/$FOLDER_ID

echo "Files sucessfully donwloaded and moved to folder $FOLDER_TO_SAVE"
