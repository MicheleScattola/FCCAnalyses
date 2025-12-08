# FCCAnalyses

[![DOI](https://zenodo.org/badge/177151745.svg)](https://zenodo.org/doi/10.5281/zenodo.4767810)

![test](https://github.com/HEP-FCC/FCCAnalyses/actions/workflows/test.yml/badge.svg?branch=pre-edm4hep1)
![docs](https://github.com/HEP-FCC/FCCAnalyses/actions/workflows/docs.yml/badge.svg?branch=pre-edm4hep1)

Common framework for FCC related analyses. This framework allows one to write
full analysis, taking [EDM4hep](https://github.com/key4hep/EDM4hep) input ROOT
files and producing the plots.
## Z -> tau tau analysis
The overall code structure is the following:
- personalized functions are stored in `FCCAnalyses/analyzer/dataframe/src/Ztautau.cc` and `FCCAnalyses/analyzer/dataframe/FCCAnalyses/Ztautau.h`
- `/FCCAnalyses/scripts/` contains the preliminary script `analysis.py` which collects necessary information for the optimal variables
- `matrix*.C` files can be exectued to plot the confusion matrix in its simple, purity or efficiency form
- `thrust.C` can be executed to plot cos(theta) and phi of the reconstructed Thrust axis
- files in `FCCAnalyses/scripts/1prong/` contain the macros for 1-prong decays analysis. One can also run `templates.py` and `templates.C` to create and store histograms to re-weight templates of pion decays
- files in `FCCAnalyses/scripts/sel1/` provide a first selection for E>2 GeV in pion decays and the associated macros to re-plot the confusion matrices
- `/FCCAnalyses/Ztautau/plots` contains all the created plots

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

## Test 1
