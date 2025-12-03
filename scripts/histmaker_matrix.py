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
outputDir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/histmaker"
#input directory
inputDir    = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/treemaker"

# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = -1
# scale the histograms with the cross-section and integrated luminosityroot- 
doScale = False
intLumi = 5000000 # 5 /ab


# define some binning for various histograms
bins = (90, 0, 45)
binst = (50, 0.7, 1)

def build_graph(df, dataset):

    results = []
    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")
     
    # baseline histograms
    results.append(df.Histo1D(("RhoAsPi_phe", "Energy_{#gamma}", *bins), "RhoAsPi_phe"))
    results.append(df.Histo1D(("RhoAsOt_phe", "Energy_{#gamma}", *bins), "RhoAsOt_phe"))
    results.append(df.Histo1D(("OtAsRho_phe", "Energy_{#gamma}", *bins), "OtAsRho_phe"))
    
    results.append(df.Histo1D(("RhoAsPi_pht", "costheta", *binst), "RhoAsPi_pht"))
    results.append(df.Histo1D(("RhoAsOt_pht", "costheta", *binst), "RhoAsOt_pht"))
    results.append(df.Histo1D(("OtAsRho_pht", "costheta", *binst), "OtAsRho_pht"))
    
    results.append(df.Histo1D(("RhoAsPi_e", "Energy_{#gamma}", *bins), "RhoAsPi_e"))
    results.append(df.Histo1D(("ElAsPi_e", "Energy_{#gamma}", *bins), "ElAsPi_e"))
    results.append(df.Histo1D(("MuAsPi_e", "Energy_{#gamma}", *bins), "MuAsPi_e"))
    results.append(df.Histo1D(("PiAsPi_e", "Energy_{#gamma}", *bins), "PiAsPi_e"))
    results.append(df.Histo1D(("OtAsPi_e", "Energy_{#gamma}", *bins), "OtAsPi_e"))
    
    results.append(df.Histo1D(("RhoAsPi_t", "costheta", *binst), "RhoAsPi_t"))
    results.append(df.Histo1D(("ElAsPi_t", "costheta", *binst), "ElAsPi_t"))
    results.append(df.Histo1D(("MuAsPi_t", "costheta", *binst), "MuAsPi_t"))
    results.append(df.Histo1D(("PiAsPi_t", "costheta", *binst), "PiAsPi_t"))
    results.append(df.Histo1D(("OtAsPi_t", "costheta", *binst), "OtAsPi_t"))

    return results, weightsum

