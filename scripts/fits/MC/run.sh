#!/bin/bash

# Stop on error
set -e 


echo ">>> 1. Running Templates ..."
fccanalysis run templatesMC.py 

echo ">>> 3. Running Data Analysis ..."
fccanalysis run analysisMC.py

echo ">>> DONE."
