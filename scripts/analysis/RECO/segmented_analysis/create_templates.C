#include "ROOT/RDataFrame.hxx"
#include "TH1D.h"
#include "TFile.h"
#include <iostream>
#include <string>

void create_and_save(ROOT::RDF::RNode df, 
                     TFile* fOut,
                     const std::string& var_sgn, 
                     const std::string& w_plus, 
                     const std::string& w_minus, 
                     const std::string& label, 
                     const std::string& suffix, 
                     int nBins, double xMin, double xMax) {

    
    // Even events -> Use for H+ template
    // Odd events  -> Use for H- template
    // This ensures ZERO statistical correlation between the two histograms.
    
    auto df_even = df.Filter("rdfentry_ % 2 == 0", "Split Sample (Even)");
    auto df_odd  = df.Filter("rdfentry_ % 2 != 0", "Split Sample (Odd)");

    // H+ filled from Even events
    auto h_plus_ptr  = df_even.Histo1D({("h_plus_"+suffix).c_str(),  ("Template Helicity +1 (" + label + ");x;Events").c_str(), nBins, xMin, xMax}, var_sgn, w_plus);
    
    // H- filled from Odd events
    auto h_minus_ptr = df_odd.Histo1D({("h_minus_"+suffix).c_str(), ("Template Helicity -1 (" + label + ");x;Events").c_str(), nBins, xMin, xMax}, var_sgn, w_minus);

    // Clone and detach from RDataFrame
    TH1D *h_plus  = (TH1D*)h_plus_ptr->Clone(("h_template_"+suffix+"_plus").c_str());
    TH1D *h_minus = (TH1D*)h_minus_ptr->Clone(("h_template_"+suffix+"_minus").c_str());

    // Print histogram entry counts
    std::cout << "[INFO] " << label << " - h_plus entries: " << h_plus->GetEntries() << ", h_minus entries: " << h_minus->GetEntries() << std::endl;

    // Safety check for empty histograms (NaN protection)
    if (h_plus->Integral() <= 0 || h_minus->Integral() <= 0) {
        std::cerr << "[WARNING] Empty histogram for " << label << ". Skipping normalization." << std::endl;
    } else {
        // Normalization
        h_plus->Scale(1.0 / h_plus->Integral());
        h_minus->Scale(1.0 / h_minus->Integral());
    }

    // Write to file
    fOut->cd();
    h_plus->Write();
    h_minus->Write();
    
    std::cout << "[INFO] " << label << " histograms saved." << std::endl;
}

void create_templates() {
        
    // Settings
    const char* infile = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/segmented_analysis/templates_merged.root";
    const char* outdir = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/segmented_analysis/";
    
    std::string treeName = "events"; 
    
    int nBins = 38;
    double xMin = 0.05;
    double xMax = 1.0; 

    // Open Dataframe
    ROOT::EnableImplicitMT(); // Enable multi-threading for speed
    ROOT::RDataFrame df(treeName, infile);

    // Open Output ROOT File
    TString rootOutName = TString(outdir) + "templates_histograms.root";
    TFile *fOut = new TFile(rootOutName, "RECREATE");

    // Processing Pions
    std::cout << "Processing Pions..." << std::endl;
    create_and_save(df, fOut, "pi_sgn", "w_plus_pi", "w_minus_pi", "Pions", "pi", nBins, xMin, xMax);

    // Processing Muons
    std::cout << "Processing Muons..." << std::endl;
    create_and_save(df, fOut, "mu_sgn", "w_plus_mu", "w_minus_mu", "Muons", "mu", nBins, xMin, xMax);

    // Processing Electrons
    std::cout << "Processing Electrons..." << std::endl;
    create_and_save(df, fOut, "el_sgn", "w_plus_el", "w_minus_el", "Electrons", "el", nBins, xMin, xMax);

    // Processing Rho 
    // IMPORTANT: Add NaN filter for Rho to prevent the "Empty Histogram" crash
    std::cout << "Processing Rho..." << std::endl;
    auto df_rho_clean = df.Filter("!std::isnan(rho_sgn[0])", "NaN Filter Rho");
    create_and_save(df_rho_clean, fOut, "rho_sgn", "w_plus_rho", "w_minus_rho", "Rho", "rho", 40, -1., 1.);

    // Close file
    fOut->Close();

    std::cout << "\n[INFO] All template histograms created and saved in: " << rootOutName << std::endl;
}
