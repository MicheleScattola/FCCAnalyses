#!/bin/bash

# Arguments passed by Condor (from the manifest file)
# $1 = Path to the text file containing the list of ROOT files
# $2 = Path where the output ROOT file should be saved
# $3 = The python script to run (templatesMC.py or analysisMC.py)

INPUT_LIST=$1
OUTPUT_FILE=$2
SCRIPT_NAME=$3

echo ">>> Job started on" `hostname`
echo "Processing list: $INPUT_LIST"
echo "Running script:  $SCRIPT_NAME"
echo "Output target:   $OUTPUT_FILE"

# Setup environment (adjust if needed)
source /afs/cern.ch/user/s/scattola/FCCAnalyses/setup.sh

# --- EXECUTION WITH XARGS ---
# We cat the list file and pipe it to xargs.
# xargs appends the file names to the end of the fccanalysis command.

cat $INPUT_LIST | xargs fccanalysis run $SCRIPT_NAME \
    --output $OUTPUT_FILE \
    --n-threads 1 \
    --files-list 

echo ">>> Job Complete"