#include "FCCAnalyses/Fitter.h"
#include <string>
#include <iostream>
#include <vector>
#include <glob.h>
#include <cmath>
#include <cstdio>
#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TColor.h>
#include <TStyle.h>
#include <TLine.h>
#include <TText.h>
#include <TLatex.h>
#include <TLegend.h>
#include <sys/stat.h>
#include <fstream>
#include <iomanip>


void fitALL() {

    // =========================================================================
    // 1. CONFIGURATION
    // =========================================================================
    const std::string merged_data_file = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/angular_analysis/output_merged.root";
    const std::string template_dir = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/angular_analysis/";
    const std::string output_dir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/RECO/";
    const std::string treeName = "events";

    // RP_costheta binning settings
    double cos_min = -1.0;
    double cos_max = 1.0;
    double step = 0.2;

    std::cout << ">>> Starting Polarization Fits on RP_costheta bins..." << std::endl;
    std::cout << ">>> Using merged data file: " << merged_data_file << std::endl;

    // =========================================================================
    // 2. STORAGE FOR RESULTS: For each bin, store single P and P_err values
    // =========================================================================
    struct BinResults {
        double P_el, P_mu, P_pi, P_rho;
        double err_el, err_mu, err_pi, err_rho;
        
        BinResults() : P_el(0), P_mu(0), P_pi(0), P_rho(0), 
                       err_el(0), err_mu(0), err_pi(0), err_rho(0) {}
    };
    
    const int n_bins = static_cast<int>((cos_max - cos_min) / step);
    std::vector<BinResults> results_per_bin(n_bins);

    // =========================================================================
    // 3. LOOP OVER RP_costheta BINS
    // =========================================================================
    
    int bin_index = 0;
    for (double bin_low = cos_min; bin_low < cos_max; bin_low += step) {
        double bin_high = bin_low + step;
        
        std::cout << "\n========================================" << std::endl;
        std::cout << "Processing RP_costheta bin [" << bin_low << ", " << bin_high << ") - Index: " << bin_index << std::endl;
        std::cout << "========================================" << std::endl;

        // Load the corresponding template file for this bin
        TString infile_templates = Form("%stemplates_histograms_%d.root", template_dir.c_str(), bin_index);
        std::cout << "  Loading template file: " << infile_templates << std::endl;
        
        BinResults& bin_results = results_per_bin[bin_index];

        std::string bin_str = std::to_string(bin_index);
        std::string el_col = "el_sgn_bin" + bin_str;
        std::string mu_col = "mu_sgn_bin" + bin_str;
        std::string pi_col = "pi_sgn_bin" + bin_str;
        std::string rho_col = "rho_sgn_bin" + bin_str;

        // ===================================================================
        // 4. FIT CHANNELS WITH FILTERED DATA
        // ===================================================================

        // --- Electron Channel ---
        std::cout << "  Fitting Electron channel..." << std::endl;
        Fitter::myFit result_el = Fitter::fit_filtered(
            merged_data_file, infile_templates.Data(),
            treeName, el_col, 
            "h_template_el_plus", 
            "h_template_el_minus"
        );
        if (result_el.success) {
            bin_results.P_el = result_el.P_tau;
            bin_results.err_el = result_el.P_err;
            std::cout << "    Electron: P = " << result_el.P_tau << " ± " << result_el.P_err << std::endl;
        } else {
            std::cout << "    Electron: Fit failed" << std::endl;
        }

        // --- Muon Channel ---
        std::cout << "  Fitting Muon channel..." << std::endl;
        Fitter::myFit result_mu = Fitter::fit_filtered(
            merged_data_file, infile_templates.Data(),
            treeName, mu_col, 
            "h_template_mu_plus", 
            "h_template_mu_minus"
        );
        if (result_mu.success) {
            bin_results.P_mu = result_mu.P_tau;
            bin_results.err_mu = result_mu.P_err;
            std::cout << "    Muon: P = " << result_mu.P_tau << " ± " << result_mu.P_err << std::endl;
        } else {
            std::cout << "    Muon: Fit failed" << std::endl;
        }

        // --- Pion Channel ---
        std::cout << "  Fitting Pion channel..." << std::endl;
        Fitter::myFit result_pi = Fitter::fit_filtered(
            merged_data_file, infile_templates.Data(),
            treeName, pi_col, 
            "h_template_pi_plus", 
            "h_template_pi_minus"
        );
        if (result_pi.success) {
            bin_results.P_pi = result_pi.P_tau;
            bin_results.err_pi = result_pi.P_err;
            std::cout << "    Pion: P = " << result_pi.P_tau << " ± " << result_pi.P_err << std::endl;
        } else {
            std::cout << "    Pion: Fit failed" << std::endl;
        }

        // --- Rho Channel ---
        std::cout << "  Fitting Rho channel..." << std::endl;
        Fitter::myFit result_rho = Fitter::fit_filtered(
            merged_data_file, infile_templates.Data(),
            treeName, rho_col, 
            "h_template_rho_plus", 
            "h_template_rho_minus"
        );
        if (result_rho.success) {
            bin_results.P_rho = result_rho.P_tau;
            bin_results.err_rho = result_rho.P_err;
            std::cout << "    Rho: P = " << result_rho.P_tau << " ± " << result_rho.P_err << std::endl;
        } else {
            std::cout << "    Rho: Fit failed" << std::endl;
        }
        
        std::cout << "  Bin " << bin_index << " processing completed." << std::endl;
        bin_index++;
    }

    // =========================================================================
    // 5. COMBINE CHANNELS AND CREATE PLOT
    // =========================================================================

    // Ensure output directory exists
    mkdir(output_dir.c_str(), 0755);

    // Prepare data for TGraphErrors
    const int n_bins_plot = static_cast<int>(results_per_bin.size());
    std::vector<double> x_values, y_values, x_errors, y_errors;

    std::cout << "\n\n========================================" << std::endl;
    std::cout << "COMBINED RESULTS BY RP_costheta BIN" << std::endl;
    std::cout << "========================================\n" << std::endl;

    
    // Optional: use stddev from segmented analysis instead of fit errors
    TFile* fAvg = TFile::Open("/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/segmented_analysis/average.root", "READ");
    TH1D* h_el = fAvg ? (TH1D*)fAvg->Get("h_Electron") : nullptr;
    TH1D* h_mu = fAvg ? (TH1D*)fAvg->Get("h_Muon") : nullptr;
    TH1D* h_pi = fAvg ? (TH1D*)fAvg->Get("h_Pion") : nullptr;
    TH1D* h_rho = fAvg ? (TH1D*)fAvg->Get("h_Rho") : nullptr;
    double std_el = h_el ? h_el->GetStdDev() : 0.0;
    double std_mu = h_mu ? h_mu->GetStdDev() : 0.0;
    double std_pi = h_pi ? h_pi->GetStdDev() : 0.0;
    double std_rho = h_rho ? h_rho->GetStdDev() : 0.0;
    if (fAvg) fAvg->Close();
    

    bin_index = 0;
    for (double bin_low = cos_min; bin_low < cos_max; bin_low += step) {
        double bin_high = bin_low + step;
        double bin_center = (bin_low + bin_high) / 2.0;
        double bin_width = (bin_high - bin_low) / 2.0; // half width for error bars

        const BinResults& bin_results = results_per_bin[bin_index];

        // Weighted average of all channels using inverse variance weighting
        // w_i = 1 / sigma_i^2
        std::vector<double> P_values = {bin_results.P_el, bin_results.P_mu, bin_results.P_pi, bin_results.P_rho};
        std::vector<double> P_errors = {bin_results.err_el, bin_results.err_mu, bin_results.err_pi, bin_results.err_rho};
        
        // Optional: override per-bin fit errors with stddev from average.root
        P_errors = {std_el, std_mu, std_pi, std_rho};
        
        
        double sum_weights = 0.0;
        double sum_weighted_P = 0.0;
        int valid_channels = 0;

        for (size_t i = 0; i < P_values.size(); ++i) {
            if (P_errors[i] > 0) {
                double weight = 1.0 / (P_errors[i] * P_errors[i]);
                sum_weighted_P += P_values[i] * weight;
                sum_weights += weight;
                valid_channels++;
            }
        }

        if (sum_weights > 0 && valid_channels > 0) {
            double P_combined = sum_weighted_P / sum_weights;
            double P_combined_err = sqrt(1.0 / sum_weights);

            x_values.push_back(bin_center);
            y_values.push_back(P_combined);
            x_errors.push_back(bin_width);
            y_errors.push_back(P_combined_err);

            std::cout << "Bin " << bin_index << " [" << bin_low << ", " << bin_high << "):" << std::endl;
            std::cout << "  P_el = " << std::fixed << std::setprecision(6) << bin_results.P_el << " ± " << bin_results.err_el << std::endl;
            std::cout << "  P_mu = " << std::fixed << std::setprecision(6) << bin_results.P_mu << " ± " << bin_results.err_mu << std::endl;
            std::cout << "  P_pi = " << std::fixed << std::setprecision(6) << bin_results.P_pi << " ± " << bin_results.err_pi << std::endl;
            std::cout << "  P_rho = " << std::fixed << std::setprecision(6) << bin_results.P_rho << " ± " << bin_results.err_rho << std::endl;
            std::cout << "  Combined P = " << std::fixed << std::setprecision(6) << P_combined 
                      << " ± " << P_combined_err << " (from " << valid_channels << " channels)" << std::endl;
        } else {
            std::cout << "Bin " << bin_index << " [" << bin_low << ", " << bin_high << "): No valid results" << std::endl;
        }

        bin_index++;
    }

    std::cout << "========================================\n" << std::endl;

    // =========================================================================
    // 6. CREATE AND SAVE ROOT PLOT
    // =========================================================================

    if (x_values.empty()) {
        std::cerr << "[ERROR] No valid results to plot!" << std::endl;
        return;
    }

    TGraphErrors* graph = new TGraphErrors(x_values.size(), 
                                           x_values.data(), 
                                           y_values.data(), 
                                           x_errors.data(), 
                                           y_errors.data());

    // Also save as ROOT file
    const std::string output_root = template_dir + "universality.root";
    TFile* fout = new TFile(output_root.c_str(), "RECREATE");
    graph->Write("graph_polarization");
    fout->Close();
    std::cout << "[INFO] ROOT file saved to: " << output_root << std::endl;

    delete graph;
    std::cout << ">>> All fits completed." << std::endl;
}
