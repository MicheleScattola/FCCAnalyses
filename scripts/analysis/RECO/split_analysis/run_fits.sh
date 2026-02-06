#!/bin/bash

# Stop on error
set -e 

echo ">>> 1. Merging template files..."

python3 merge.py

echo ">>> 2. Creating template histograms..."

root -l -b -q templates.C

echo ">>> 3. Fitting all data files..."

root -l -b -q fitALL.C

echo ">>> 4. Generating confusion matrices..."

root -l -b -q ../confusion/matrix.C
root -l -b -q ../confusion/matrix_eff.C
root -l -b -q ../confusion/matrix_pur.C

echo ">>> 5. Generating invariant mass plots..."

root -l -b -q ../additional/mass.C

#echo ">>> 6. calculating rejecetion efficiencies.."

root -l -b -q rho_rejection_efficiency.C


echo ">>> DONE."
