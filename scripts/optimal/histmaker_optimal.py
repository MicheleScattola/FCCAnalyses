import math
#List of processes
processList = {
    'p8_ee_Ztautau_ecm91':{}
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
#prodTag     = "FCCee/winter2023/IDEA/"
# Link to the dictonary that contains all the cross section informations etc... (mandatory)
procDict = "FCCee_procDict_winter2023_IDEA.json"

#output directory
outputDir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/histmaker/optimal"
#input directory
inputDir    = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/treemaker/optimal"

# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = -1
# scale the histograms with the cross-section and integrated luminosityroot- 
doScale = False
intLumi = 5000000 # 5 /ab


# define some binning for various histograms
bins = (38, 0.05, 1)
bins2 = (46, 0.05, 1.2)
binst = (50, 0.7, 1)

def build_graph(df, dataset):

    results = []
    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")
     
    # baseline histograms
    results.append(df.Histo1D(("mu_sgn", "x_mu", *bins), "mu_sgn"))
    results.append(df.Histo1D(("el_sgn", "x_el", *bins), "el_sgn"))
    
    results.append(df.Histo1D(("pi_bkg", "pi_bkg", *bins), "pi_bkg"))
    results.append(df.Histo1D(("pi_sgn", "x_el", *bins), "pi_sgn"))

    return results, weightsum

