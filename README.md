
# $Z \rightarrow \tau \tau$ analysis
Analysis performed by Michele Scattola, if you have any questions don't hesitate to write at michele.scattola@studenti.unimi.it

The overall code structure is the following:
- personalized functions are stored in `FCCAnalyses/analyzer/dataframe/src/Ztautau.cc` and `FCCAnalyses/analyzer/dataframe/FCCAnalyses/Ztautau.h` . Further specifications are explained in the file **`FUNCTIONS.md`** .
- fitting functions are stored in `FCCAnalyses/analyzer/dataframe/src/Fitter.cc` and `FCCAnalyses/analyzer/dataframe/FCCAnalyses/Fitter.h` .
The methods applied are explained in the file **`FITTING.md`** .
- analysis scripts are stored in `FCCAnalyses/scripts/analysis` in folders `/MC` and `/RECO` for Monte Carlo or Reconstructed analysis.
- both analysis folder feature a `/split_analysis` section where the templates and data events are collected from indipendent event samples.
- `/FCCAnalyses/Ztautau/plots` contains all the created plots for MC and RECO.

### Usage:
**Tip:** in order to understand more easily the workflow I suggest starting from the bash scripts `run_fits.sh` which are present in most of the analysis folders.

As with most of the code it won't easily be applicable outside of the FCCAnalyses setup. The code was developed for a specific *key4hep* stack version, which should always be the version sourced on lxplus machines. The version depends on the dataset used, but it should be clearly labelled:
https://fcc-physics-events.web.cern.ch/fcc-ee/rec/winter2023/IDEA
For my `p8_ee_Ztautau_ecm91` data the stack is obtainable via
```bash
source /cvmfs/sw.hsf.org/spackages6/key4hep-stack/2022-12-23/x86_64-centos7-gcc11.2.0-opt/ll3gi/setup.sh
```

*JUNE 2026 NOTE:*
It seems that the update to AlmaLinux9 machines broke down the software stack used prior.
In order to compile now clean up the local install open up a **new** shell and source an intermediate version:
```bash
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-03-10
```
cleanup the old install and re-build
```bash
cd ~/path/to/FCCAnalyses
rm -rf build/
rm -rf install/
mkdir build install
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 4
```

I highly suggest to go through the tutorial for a better understanding of the FCCAnalyses framework ([FCC tutorial](https://hep-fcc.github.io/fcc-tutorials/main/index.html)).



# RECO analysis
## /split_analysis:
This folder produces an analysis over $20 \times 10^6$ events for both templates and data. The events are taken from indipendent samples via python scripts and convenient condor queuing system.
## /segmented_analysis:
This folder produces a repeated analysis with $20 \times 10^6$ events for templates and $1 \times 10^6$ events for data. The polarization analysis is perfomed over 80 different sets of data, in order to calculate the standard deviation of the fitting results.
## /angular_analysis:
This folder produces an analysis of the datasets binned by angular production, providing the fits for the leptonic universality test.
## /confusion:
This folder produces the plots for confusion matrices and background contamination of the reconstruction performed.

# MC analysis
This section provides the same analysis ad 'RECO/split_analysis/' but on Montecarlo data.

======================== end of my work

# FCCAnalyses

[![DOI](https://zenodo.org/badge/177151745.svg)](https://zenodo.org/doi/10.5281/zenodo.4767810)

![test](https://github.com/HEP-FCC/FCCAnalyses/actions/workflows/test.yml/badge.svg?branch=pre-edm4hep1)
![docs](https://github.com/HEP-FCC/FCCAnalyses/actions/workflows/docs.yml/badge.svg?branch=pre-edm4hep1)

Common framework for FCC related analyses. This framework allows one to write
full analysis, taking [EDM4hep](https://github.com/key4hep/EDM4hep) input ROOT
files and producing the plots.

## Quick start

In order to run over pre-generated samples from `winter2023` or `spring2021`
campaigns one needs to compile `pre-edm4hep1` branch of the FCCAnalyses in the
`2024-03-10` release (hard coded into the `setup.sh` script)

```sh
git clone --branch pre-edm4hep1 git@github.com:HEP-FCC/FCCAnalyses.git
cd FCCAnalyses
source ./setup.sh
fccanalysis build -j 8
```

To have access to the FCC pre-generated samples, one needs to be subscribed
the `fcc-eos-access` e-group (with owner approval).

Detailed documentation can be found at the [FCCAnalyses](https://hep-fcc.github.io/FCCAnalyses/) webpage.

All sample information, including Key4hep stack used for the campaign, is collected at the
[FCC Physics Events](http://fcc-physics-events.web.cern.ch/fcc-physics-events/) website.


## Contributing

As usual, if you aim at contributing to the repository, please fork it, develop your feature/analysis and submit a pull requests.

### Code formating

The preferred style of the C++ code in the FCCAnalyses is LLVM which is checked
by CI job.

To apply formatting to a given file:
```
clang-format -i -style=file /path/to/file.cpp
```
