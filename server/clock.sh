#!/bin/bash

# Infinite loop to run the squeue command every 5 seconds
while true
do
    echo "Checking squeue"
    squeue -u elopatuhina

    # Wait
    sleep 5s
done
