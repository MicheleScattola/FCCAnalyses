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
bins2 = (60, 35, 55)
bins3 = (220, 0, 220)

def build_graph(df, dataset):

    results = []
    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")
     
    # baseline histograms
    results.append(df.Histo1D(("gammagamma_e", "Energy_{#gamma}", *bins), "gammagamma_e"))
    results.append(df.Histo1D(("photons_e", "Energy_{#gamma}", *bins), "photons_e"))
    results.append(df.Histo1D(("MC_gammagamma_e", "Energy_{#gamma}", *bins), "MC_gammagamma_e"))
    results.append(df.Histo1D(("pi1gamma1_true", "Energy_{#gamma}", *bins), "pi1gamma1_true"))
    results.append(df.Histo1D(("pi1gamma1_false", "Energy_{#gamma}", *bins), "pi1gamma1_false"))
    results.append(df.Histo1D(("pi1gamma2_true", "Energy_{#gamma}", *bins), "pi1gamma2_true"))
    results.append(df.Histo1D(("pi1gamma2_false", "Energy_{#gamma}", *bins), "pi1gamma2_false"))
    results.append(df.Histo1D(("wrongrho_1ph", "Energy_{#gamma}", *bins), "wrongrho_1ph"))
    results.append(df.Histo1D(("wrongrho_2ph", "Energy_{#gamma}", *bins), "wrongrho_2ph"))
    results.append(df.Histo1D(("sum_wrong2ph", "Energy_{#gamma}", *bins), "sum_wrong2ph"))
    results.append(df.Histo1D(("sum_true2ph", "Energy_{#gamma}", *bins), "sum_true2ph"))
    results.append(df.Histo1D(("taus_e_collinear", "Energy_{#tau}", *bins2), "taus_e_collinear"))
    results.append(df.Histo1D(("parentPDG1", "PDG1", *bins3), "parentPDG1"))
    results.append(df.Histo1D(("parentPDG2", "PDG2", *bins3), "parentPDG2"))
    
    results.append(df.Histo1D(("el1g_all", "Energy_{#gamma}", *bins), "el1g_all"))
    results.append(df.Histo1D(("mu1g_all", "Energy_{#gamma}", *bins), "mu1g_all"))
    results.append(df.Histo1D(("tau1g_all", "Energy_{#gamma}", *bins), "tau1g_all"))
    results.append(df.Histo1D(("ph1g_all", "Energy_{#gamma}", *bins), "ph1g_all"))
    results.append(df.Histo1D(("pi01g_all", "Energy_{#gamma}", *bins), "pi01g_all"))
    results.append(df.Histo1D(("el1g_true", "Energy_{#gamma}", *bins), "el1g_true"))
    results.append(df.Histo1D(("mu1g_true", "Energy_{#gamma}", *bins), "mu1g_true"))
    results.append(df.Histo1D(("tau1g_true", "Energy_{#gamma}", *bins), "tau1g_true"))
    results.append(df.Histo1D(("ph1g_true", "Energy_{#gamma}", *bins), "ph1g_true"))
    results.append(df.Histo1D(("pi01g_true", "Energy_{#gamma}", *bins), "pi01g_true"))
    

    return results, weightsum

