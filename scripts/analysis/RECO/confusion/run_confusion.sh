#!/bin/bash

# Stop on error
set -e 
echo ">>> 1. Generating confusion matrices..."

root -l -b -q confusion/matrix.C
root -l -b -q confusion/matrix_eff.C
root -l -b -q confusion/matrix_pur.C

echo ">>> DONE."