#!/bin/bash

# Stop on error
set -e 


echo ">>> 1. Running Templates ... submitting jobs to batch system ..."
fccanalysis run templatesMC.py 

echo ">>> 2. Running Data Analysis ... submitting jobs to batch system ..."
fccanalysis run analysisMC.py

echo ">>> DONE."
