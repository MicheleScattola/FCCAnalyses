#!/bin/bash

# Stop on error
set -e 

echo ">>> 1. Generating templates..."

root -l -b -q templates.C

echo ">>> 2. Fitting all decay channels..."

root -l -b -q fitALL.C

echo ">>> DONE."