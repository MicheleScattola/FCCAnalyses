#!/bin/bash

# Stop on error
set -e 


echo ">>> 1. Running Templates ... submitting jobs to batch system ..."
fccanalysis run templates.py 

echo ">>> 2. Running Data Analysis ... submitting jobs to batch system ..."
fccanalysis run analysis.py

echo ">>> DONE."
