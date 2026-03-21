#!/bin/bash

for val in {2..13}; do
    echo "====================================="
    echo "🚀 RUNNING SIMULATION NO. : $val"
    echo "📁 Using script: run_all.sh"
    echo "🕒 Started at: $(date)"
    echo "====================================="

    # Send the number, then an extra newline for the PAUSE
    (echo $val; echo) | ./a.out

    echo "----------------------"
done

