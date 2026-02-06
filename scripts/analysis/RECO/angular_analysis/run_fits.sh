#!/bin/bash

# Stop on error
set -e 

echo ">>> 1. Merging template files..."

python3 merge.py

echo ">>> 2. Creating template histograms..."

root -l -b -q create_templates.C

echo ">>> 3. Fitting all data files..."

root -l -b -q fitALL.C

echo ">>> 4. Creating summary plots..."

root -l -b -q plot.C

echo ">>> DONE."
