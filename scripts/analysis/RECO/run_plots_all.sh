#!/bin/bash

# Stop on error
set -e 

echo ">>> 1. Generating templates..."

root -l -b -q templates.C

echo ">>> 2. Fitting all decay channels..."

root -l -b -q fitALL.C

echo ">>> 3. Generating confusion matrices..."

root -l -b -q confusion/matrix.C
root -l -b -q confusion/matrix_eff.C
root -l -b -q confusion/matrix_pur.C

echo ">>> 4. Generating additional plots..."

root -l -b -q additional/thrust.C
root -l -b -q additional/mass.C

echo ">>> DONE."