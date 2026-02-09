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
    const std::string output_dir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/RECO/angular_analysis/";
    const std::string treeName = "events";

    // RP_costheta binning settings
    const double cos_min = -0.95;
    const double cos_max = 0.95;
    const int n_costheta_bins = 20;
    const double step = (cos_max - cos_min) / n_costheta_bins;

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
    
    std::vector<BinResults> results_per_bin(n_costheta_bins);

    // =========================================================================
    // 3. LOOP OVER RP_costheta BINS
    // =========================================================================
    
    for (int bin_index = 0; bin_index < n_costheta_bins; ++bin_index) {
        double bin_low = cos_min + bin_index * step;
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
            "h_template_el_minus",
            true, output_dir, bin_low, bin_high
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
            "h_template_mu_minus",
            true, output_dir, bin_low, bin_high
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
            "h_template_pi_minus",
            true, output_dir, bin_low, bin_high
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
            "h_template_rho_minus",
            true, output_dir, bin_low, bin_high
        );
        if (result_rho.success) {
            bin_results.P_rho = result_rho.P_tau;
            bin_results.err_rho = result_rho.P_err;
            std::cout << "    Rho: P = " << result_rho.P_tau << " ± " << result_rho.P_err << std::endl;
        } else {
            std::cout << "    Rho: Fit failed" << std::endl;
        }
        
        std::cout << "  Bin " << bin_index << " processing completed." << std::endl;
    }

    // =========================================================================
    // 5. COMBINE CHANNELS AND CREATE PLOT
    // =========================================================================

    // Ensure output directory exists
    mkdir(output_dir.c_str(), 0755);

    // Prepare data for TGraphErrors
    const int n_bins_plot = static_cast<int>(results_per_bin.size());
    std::vector<double> x_values, y_values, x_errors, y_errors;
    
    // Separate vectors for each channel
    std::vector<double> x_values_el, y_values_el, x_errors_el, y_errors_el;
    std::vector<double> x_values_mu, y_values_mu, x_errors_mu, y_errors_mu;
    std::vector<double> x_values_pi, y_values_pi, x_errors_pi, y_errors_pi;
    std::vector<double> x_values_rho, y_values_rho, x_errors_rho, y_errors_rho;

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
    const bool use_segmented_errors = (std_el > 0 && std_mu > 0 && std_pi > 0 && std_rho > 0);
    if (fAvg) fAvg->Close();
    

    for (int bin_index = 0; bin_index < n_costheta_bins; ++bin_index) {
        double bin_low = cos_min + bin_index * step;
        double bin_high = bin_low + step;
        double bin_center = (bin_low + bin_high) / 2.0;
        double bin_width = (bin_high - bin_low) / 2.0; // half width for error bars

        const BinResults& bin_results = results_per_bin[bin_index];

        // Weighted average of all channels using inverse variance weighting
        // w_i = 1 / sigma_i^2
        std::vector<double> P_values = {bin_results.P_el, bin_results.P_mu, bin_results.P_pi, bin_results.P_rho};
        std::vector<double> P_errors = {bin_results.err_el, bin_results.err_mu, bin_results.err_pi, bin_results.err_rho};
        
        // Use segmented analysis errors for the quadratic sum in the weighted average
        /*if (use_segmented_errors) {
            P_errors = {std_el, std_mu, std_pi, std_rho};
        }*/
        
        // Store individual channel results
        if (bin_results.err_el > 0) {
            x_values_el.push_back(bin_center);
            y_values_el.push_back(bin_results.P_el);
            x_errors_el.push_back(bin_width);
            y_errors_el.push_back(bin_results.err_el);
        }
        if (bin_results.err_mu > 0) {
            x_values_mu.push_back(bin_center);
            y_values_mu.push_back(bin_results.P_mu);
            x_errors_mu.push_back(bin_width);
            y_errors_mu.push_back(bin_results.err_mu);
        }
        if (bin_results.err_pi > 0) {
            x_values_pi.push_back(bin_center);
            y_values_pi.push_back(bin_results.P_pi);
            x_errors_pi.push_back(bin_width);
            y_errors_pi.push_back(bin_results.err_pi);
        }
        if (bin_results.err_rho > 0) {
            x_values_rho.push_back(bin_center);
            y_values_rho.push_back(bin_results.P_rho);
            x_errors_rho.push_back(bin_width);
            y_errors_rho.push_back(bin_results.err_rho);
        }
        
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
                      << " ± " << P_combined_err << " (from " << valid_channels << " channels, "
                      << (use_segmented_errors ? "segmented errors" : "per-bin fit errors") << ")" << std::endl;
        } else {
            std::cout << "Bin " << bin_index << " [" << bin_low << ", " << bin_high << "): No valid results" << std::endl;
        }

    }

    std::cout << "========================================\n" << std::endl;

    // =========================================================================
    // 6. CREATE AND SAVE ROOT HISTOGRAMS (individual channels + combined)
    // =========================================================================

    if (x_values.empty()) {
        std::cerr << "[ERROR] No valid results to plot!" << std::endl;
        return;
    }

    // Create TGraphErrors for each channel
    TGraphErrors* graph_el = new TGraphErrors(x_values_el.size(), 
                                               x_values_el.data(), 
                                               y_values_el.data(), 
                                               x_errors_el.data(), 
                                               y_errors_el.data());
    graph_el->SetName("el_theta");
    graph_el->SetTitle("Electron Polarization vs #cos(#theta)");

    TGraphErrors* graph_mu = new TGraphErrors(x_values_mu.size(), 
                                               x_values_mu.data(), 
                                               y_values_mu.data(), 
                                               x_errors_mu.data(), 
                                               y_errors_mu.data());
    graph_mu->SetName("mu_theta");
    graph_mu->SetTitle("Muon Polarization vs #cos(#theta)");

    TGraphErrors* graph_pi = new TGraphErrors(x_values_pi.size(), 
                                               x_values_pi.data(), 
                                               y_values_pi.data(), 
                                               x_errors_pi.data(), 
                                               y_errors_pi.data());
    graph_pi->SetName("pi_theta");
    graph_pi->SetTitle("Pion Polarization vs #cos(#theta)");

    TGraphErrors* graph_rho = new TGraphErrors(x_values_rho.size(), 
                                                x_values_rho.data(), 
                                                y_values_rho.data(), 
                                                x_errors_rho.data(), 
                                                y_errors_rho.data());
    graph_rho->SetName("rho_theta");
    graph_rho->SetTitle("Rho Polarization vs #cos(#theta)");

    // Combined graph
    TGraphErrors* graph_combined = new TGraphErrors(x_values.size(), 
                                                     x_values.data(), 
                                                     y_values.data(), 
                                                     x_errors.data(), 
                                                     y_errors.data());
    graph_combined->SetName("combined_theta");
    graph_combined->SetTitle("Combined Polarization vs #cos(#theta)");

    // Save as ROOT file
    const std::string output_root = template_dir + "universality.root";
    TFile* fout = new TFile(output_root.c_str(), "RECREATE");
    graph_el->Write();
    graph_mu->Write();
    graph_pi->Write();
    graph_rho->Write();
    graph_combined->Write();
    fout->Close();
    std::cout << "[INFO] ROOT file saved to: " << output_root << std::endl;

    delete graph_el;
    delete graph_mu;
    delete graph_pi;
    delete graph_rho;
    delete graph_combined;
    std::cout << ">>> All fits completed." << std::endl;
}
