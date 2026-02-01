#!/bin/bash

# Stop on error
set -e

echo ">>> Formatting headers in source files..."

clang-format -i -style=file /afs/cern.ch/user/s/scattola/FCCAnalyses/analyzers/dataframe/src/Ztautau.cc
clang-format -i -style=file /afs/cern.ch/user/s/scattola/FCCAnalyses/analyzers/dataframe/src/Fitter.cc

clang-format -i -style=file /afs/cern.ch/user/s/scattola/FCCAnalyses/analyzers/dataframe/FCCAnalyses/Ztautau.h
clang-format -i -style=file /afs/cern.ch/user/s/scattola/FCCAnalyses/analyzers/dataframe/FCCAnalyses/Fitter.h

echo ">>> DONE."