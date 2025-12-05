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
outputDir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/histmaker"
#input directory
inputDir    = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/treemaker"

# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = -1
# scale the histograms with the cross-section and integrated luminosityroot- 
doScale = False
intLumi = 5000000 # 5 /ab


# define some binning for various histograms
bins = (90, 0, 45)
binst = (20, 0.8, 1)

def build_graph(df, dataset):

    results = []
    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")
     
    # baseline histograms
    results.append(df.Histo1D(("RhoAsPi_e", "Energy_{#gamma}", *bins), "RhoAsPi_e"))
    results.append(df.Histo1D(("ElAsPi_e", "Energy_{#gamma}", *bins), "ElAsPi_e"))
    results.append(df.Histo1D(("MuAsPi_e", "Energy_{#gamma}", *bins), "MuAsPi_e"))
    results.append(df.Histo1D(("PiAsPi_e", "Energy_{#gamma}", *bins), "PiAsPi_e"))
    results.append(df.Histo1D(("OtAsPi_e", "Energy_{#gamma}", *bins), "OtAsPi_e"))
    results.append(df.Histo1D(("A1AsPi_e", "Energy_{#gamma}", *bins), "A1AsPi_e"))
    
    #results.append(df.Histo1D(("RhoAsPi_t", "costheta", *binst), "RhoAsPi_t"))
    #results.append(df.Histo1D(("ElAsPi_t", "costheta", *binst), "ElAsPi_t"))
    #results.append(df.Histo1D(("MuAsPi_t", "costheta", *binst), "MuAsPi_t"))
    #results.append(df.Histo1D(("PiAsPi_t", "costheta", *binst), "PiAsPi_t"))
    #results.append(df.Histo1D(("OtAsPi_t", "costheta", *binst), "OtAsPi_t"))
    #results.append(df.Histo1D(("A1AsPi_t", "costheta", *binst), "A1AsPi_t"))

    return results, weightsum

