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
    std::string size_filter = var_sgn + ".size() > 0 && " + w_plus + ".size() == " + var_sgn + ".size() && " + w_minus + ".size() == " + var_sgn + ".size()";
    auto df_valid = df.Filter(size_filter, "non-empty bin filter");
    // auto df_even = df_valid.Filter("rdfentry_ % 2 == 0", "Split Sample (Even)");
    // auto df_odd  = df_valid.Filter("rdfentry_ % 2 != 0", "Split Sample (Odd)");

    // H+ filled from all events (even and odd combined for better statistics)
    auto h_plus_ptr  = df_valid.Histo1D({("h_plus_"+suffix).c_str(),  ("Template Helicity +1 (" + label + ");x;Events").c_str(), nBins, xMin, xMax}, var_sgn, w_plus);
    
    // H- filled from all events (even and odd combined for better statistics)
    auto h_minus_ptr = df_valid.Histo1D({("h_minus_"+suffix).c_str(), ("Template Helicity -1 (" + label + ");x;Events").c_str(), nBins, xMin, xMax}, var_sgn, w_minus);

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
    const char* infile = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/angular_analysis/templates_merged.root";
    const char* outdir = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/angular_analysis/";
    
    std::string treeName = "events"; 
    
    int nBins = 38;
    double xMin = 0.05;
    double xMax = 1.0; 

    // RP_costheta binning settings
    const double cos_min = -0.95;
    const double cos_max = 0.95;
    const int n_costheta_bins = 20;
    const double step = (cos_max - cos_min) / n_costheta_bins;
    
    // Open Dataframe
    ROOT::EnableImplicitMT(); // Enable multi-threading for speed
    ROOT::RDataFrame df(treeName, infile);

    // Loop over RP_costheta bins
    for (int bin_index = 0; bin_index < n_costheta_bins; ++bin_index) {
        double bin_low = cos_min + bin_index * step;
        double bin_high = bin_low + step;
        
        std::cout << "\n========================================" << std::endl;
        std::cout << "Processing RP_costheta bin [" << bin_low << ", " << bin_high << ")" << std::endl;
        std::cout << "========================================\n" << std::endl;
        
        // Open Output ROOT File for this bin
        TString rootOutName = Form("%stemplates_histograms_%d.root", outdir, bin_index);
        TFile *fOut = new TFile(rootOutName, "RECREATE");

        std::string bin_str = std::to_string(bin_index);
        std::string pi_sgn_bin = "pi_sgn_bin" + bin_str;
        std::string mu_sgn_bin = "mu_sgn_bin" + bin_str;
        std::string el_sgn_bin = "el_sgn_bin" + bin_str;
        std::string rho_sgn_bin = "rho_sgn_bin" + bin_str;

        std::string w_plus_pi_bin = "w_plus_pi_bin" + bin_str;
        std::string w_minus_pi_bin = "w_minus_pi_bin" + bin_str;
        std::string w_plus_mu_bin = "w_plus_mu_bin" + bin_str;
        std::string w_minus_mu_bin = "w_minus_mu_bin" + bin_str;
        std::string w_plus_el_bin = "w_plus_el_bin" + bin_str;
        std::string w_minus_el_bin = "w_minus_el_bin" + bin_str;
        std::string w_plus_rho_bin = "w_plus_rho_bin" + bin_str;
        std::string w_minus_rho_bin = "w_minus_rho_bin" + bin_str;

        // Processing Pions
        std::cout << "Processing Pions..." << std::endl;
        create_and_save(df, fOut, pi_sgn_bin, w_plus_pi_bin, w_minus_pi_bin, "Pions", "pi", nBins, xMin, xMax);

        // Processing Muons
        std::cout << "Processing Muons..." << std::endl;
        create_and_save(df, fOut, mu_sgn_bin, w_plus_mu_bin, w_minus_mu_bin, "Muons", "mu", nBins, xMin, xMax);

        // Processing Electrons
        std::cout << "Processing Electrons..." << std::endl;
        create_and_save(df, fOut, el_sgn_bin, w_plus_el_bin, w_minus_el_bin, "Electrons", "el", nBins, xMin, xMax);

        // Processing Rho 
        std::cout << "Processing Rho..." << std::endl;
        create_and_save(df, fOut, rho_sgn_bin, w_plus_rho_bin, w_minus_rho_bin, "Rho", "rho", 40, -1., 1.);

        // Close file
        fOut->Close();

        std::cout << "\n[INFO] Template histograms for bin " << bin_index << " saved in: " << rootOutName << std::endl;
        
    }
    
    std::cout << "\n[INFO] All template histograms created across " << n_costheta_bins
              << " RP_costheta bins." << std::endl;
}
