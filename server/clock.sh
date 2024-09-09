#!/bin/bash

# Help function
show_help() {
        echo "Usage: $0 [--help HELP] [--time TIME] [--user USER]"
        echo
        echo  "Check squeue."
        echo
        echo "optional arguments:"
        echo
        echo " --help, -h               Show this message and exit."
        echo " --time TIME, -t TIME     Delay interval in seconds. Default value is '5s'."
        echo " --user USER, -u USER     Name of user. Default value is 'elopatuhina'."
        echo
        echo " To quite program use 'Ctrl+C'."
        echo 
}

# Set the interval for sleep (default is 5 seconds)
time="5s"
user="elopatuhina"

# Read command-line arguments
while [[ $# -gt 0 ]]; do
        case "$1" in
                -h|--help)
                        show_help
                        exit 0
                        ;;

                -t|--time)
                        if [[ -n $2 ]]; then
                                time="$2"
                        else
                                echo "Error: --time requires a non-empty argument. For details see --help."
                                exit 1
                        fi
                        ;;

                -u|--user)
                        if [[ -n $2 ]]; then
                                user="$2"
                        else
                                echo "Error: --user requires a non-empty argument. For details see --help."
                                exit 1
                        fi
                        ;;
                
                
                *)
                        echo "Invalide flag $1. For details see --help."
                        exit 1
                        ;;
        esac
        shift 2
done

# Infinite loop to run the squeue command every 5 seconds
while true
do
    echo "Checking squeue:"
    squeue -u "$user"

    # Wait
    sleep "$time"
done
