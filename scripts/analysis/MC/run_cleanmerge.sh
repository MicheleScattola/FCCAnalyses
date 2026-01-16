#!/bin/bash

# Stop on error
set -e 


echo ">>> 1. Merging templates on eos/user/s/scattola ..."
python merge_templates.py

echo ">>> 2. Cleaning template chunks ..."
rm -r /eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/MC/templates/p8_ee_Ztautau_ecm91

echo ">>> 3. Merging data on eos/user/s/scattola ..."
python merge_data.py

echo ">>> 4. Cleaning data chunks ..."
rm -r /eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/MC/p8_ee_Ztautau_ecm91

echo ">>> DONE."