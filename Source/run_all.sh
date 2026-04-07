#!/bin/bash

# Set to 1 (Sim) or 2 (DB)
FILE_TYPE=1

for val in {1..1}; do
    echo "====================================="
    echo " RUNNING SIMULATION NO. : $val"
    echo " Using script: run_all.sh"
    echo " Started at: $(date)"
    echo "====================================="
    # Send: file type choice, then sim number, then extra newline for PAUSE
    (echo $FILE_TYPE; echo $val; echo) | ./a.out
    echo "----------------------"
done
