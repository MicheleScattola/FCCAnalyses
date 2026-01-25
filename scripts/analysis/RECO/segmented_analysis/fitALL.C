#include "FCCAnalyses/Fitter.h"
#include <string>
#include <iostream>
#include <vector>
#include <glob.h>
#include <cmath>
#include <TCanvas.h>
#include <TH1D.h>
#include <TColor.h>
#include <TStyle.h>
#include <TLine.h>
#include <TText.h>

// Helper function to get all data_*.root files
std::vector<std::string> getDataFiles(const std::string& directory) {
    std::vector<std::string> files;
    glob_t glob_result;
    std::string pattern = directory + "/output_*.root";
    
    glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    
    for(unsigned int i = 0; i < glob_result.gl_pathc; ++i) {
        files.push_back(std::string(glob_result.gl_pathv[i]));
    }
    
    globfree(&glob_result);
    return files;
}

// Helper function to calculate mean and standard deviation
void calculateStats(const std::vector<double>& values, double& mean, double& stddev) {
    if (values.empty()) {
        mean = 0.0;
        stddev = 0.0;
        return;
    }
    
    double sum = 0.0;
    for (double val : values) {
        sum += val;
    }
    mean = sum / values.size();
    
    double variance = 0.0;
    for (double val : values) {
        variance += (val - mean) * (val - mean);
    }
    stddev = sqrt(variance / values.size());
}

// Build a ROOT histogram for a channel
TH1D* makeHist(const std::vector<double>& values, const std::string& name, const std::string& title, double mean, double stddev) {
    // Define range centered on mean with width of 2 sigma per side (4 sigma total)
    const double xmin = mean - 3.0 * stddev;
    const double xmax = mean + 3.0 * stddev;
    const int nbins = 30;

    TH1D* h = new TH1D(name.c_str(), title.c_str(), nbins, xmin, xmax);
    h->SetLineWidth(2);
    h->SetLineColor(kRed + 1);
    for (double v : values) {
        h->Fill(v);
    }
    h->SetStats(true);
    return h;
}

void fitALL() {

    // =========================================================================
    // 1. CONFIGURATION
    // =========================================================================
    const std::string data_dir = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/segmented_analysis/";
    const std::string output_dir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/RECO";
    const std::string infile_templates = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/segmented_analysis/templates_histograms.root";
    const std::string treeName = "events";

    std::cout << ">>> Starting Polarization Fits on Multiple Files..." << std::endl;

    // Get all data files
    std::vector<std::string> dataFiles = getDataFiles(data_dir);
    std::cout << "Found " << dataFiles.size() << " data files." << std::endl;

    if (dataFiles.empty()) {
        std::cerr << "No data files found!" << std::endl;
        return;
    }

    // =========================================================================
    // 2. STORAGE FOR RESULTS
    // =========================================================================
    std::vector<double> P_el, P_mu, P_pi, P_rho;
    std::vector<double> err_el, err_mu, err_pi, err_rho;

    // =========================================================================
    // 3. LOOP OVER DATA FILES AND FIT EACH CHANNEL
    // =========================================================================
    
    for (size_t i = 0; i < dataFiles.size(); ++i) {
        const std::string& dataFile = dataFiles[i];
        std::cout << "\n>>> Processing file " << (i+1) << "/" << dataFiles.size() << ": " << dataFile << std::endl;

        // --- Electron Channel ---
        Fitter::myFit result_el = Fitter::fit_no_plot(
            dataFile, infile_templates,
            treeName, "el_sgn", 
            "h_template_el_plus", 
            "h_template_el_minus"
        );
        if (result_el.success) {
            P_el.push_back(result_el.P_tau);
            err_el.push_back(result_el.P_err);
        }

        // --- Muon Channel ---
        Fitter::myFit result_mu = Fitter::fit_no_plot(
            dataFile, infile_templates,
            treeName, "mu_sgn", 
            "h_template_mu_plus", 
            "h_template_mu_minus"
        );
        if (result_mu.success) {
            P_mu.push_back(result_mu.P_tau);
            err_mu.push_back(result_mu.P_err);
        }

        // --- Pion Channel ---
        Fitter::myFit result_pi = Fitter::fit_no_plot(
            dataFile, infile_templates,
            treeName, "pi_sgn", 
            "h_template_pi_plus", 
            "h_template_pi_minus"
        );
        if (result_pi.success) {
            P_pi.push_back(result_pi.P_tau);
            err_pi.push_back(result_pi.P_err);
        }

        // --- Rho Channel ---
        Fitter::myFit result_rho = Fitter::fit_no_plot(
            dataFile, infile_templates,
            treeName, "rho_sgn", 
            "h_template_rho_plus", 
            "h_template_rho_minus"
        );
        if (result_rho.success) {
            P_rho.push_back(result_rho.P_tau);
            err_rho.push_back(result_rho.P_err);
        }
    }

    // =========================================================================
    // 4. CALCULATE AND DISPLAY STATISTICS
    // =========================================================================
    
    std::cout << "\n\n========================================" << std::endl;
    std::cout << "RESULTS SUMMARY" << std::endl;
    std::cout << "========================================\n" << std::endl;

    double mean_el, stddev_el, mean_mu, stddev_mu, mean_pi, stddev_pi, mean_rho, stddev_rho;

    // Electron Channel
    calculateStats(P_el, mean_el, stddev_el);
    std::cout << "Electron Channel:" << std::endl;
    std::cout << "  Number of successful fits: " << P_el.size() << std::endl;
    std::cout << "  Mean P_tau: " << mean_el << std::endl;
    std::cout << "  Std Dev:    " << stddev_el << std::endl;
    std::cout << std::endl;

    // Muon Channel
    calculateStats(P_mu, mean_mu, stddev_mu);
    std::cout << "Muon Channel:" << std::endl;
    std::cout << "  Number of successful fits: " << P_mu.size() << std::endl;
    std::cout << "  Mean P_tau: " << mean_mu << std::endl;
    std::cout << "  Std Dev:    " << stddev_mu << std::endl;
    std::cout << std::endl;

    // Pion Channel
    calculateStats(P_pi, mean_pi, stddev_pi);
    std::cout << "Pion Channel:" << std::endl;
    std::cout << "  Number of successful fits: " << P_pi.size() << std::endl;
    std::cout << "  Mean P_tau: " << mean_pi << std::endl;
    std::cout << "  Std Dev:    " << stddev_pi << std::endl;
    std::cout << std::endl;

    // Rho Channel
    calculateStats(P_rho, mean_rho, stddev_rho);
    std::cout << "Rho Channel:" << std::endl;
    std::cout << "  Number of successful fits: " << P_rho.size() << std::endl;
    std::cout << "  Mean P_tau: " << mean_rho << std::endl;
    std::cout << "  Std Dev:    " << stddev_rho << std::endl;
    std::cout << std::endl;

    std::cout << "========================================" << std::endl;

    // =========================================================================
    // 5. DRAW DISTRIBUTIONS (2x2 CANVAS) AND SAVE
    // =========================================================================

    gStyle->SetOptStat(1111); // Show entries, mean, and RMS on stat box

    TCanvas* c = new TCanvas("c", "Polarization distributions", 1200, 900);
    c->Divide(2, 2);

    std::vector<std::pair<std::vector<double>*, std::pair<double, double>>> channels = {
        {&P_el, {mean_el, stddev_el}},
        {&P_mu, {mean_mu, stddev_mu}},
        {&P_pi, {mean_pi, stddev_pi}},
        {&P_rho, {mean_rho, stddev_rho}}
    };

    std::vector<std::string> labels = {"Electron", "Muon", "Pion", "Rho"};
    const double P_tau_value = -0.14719; // Reference polarization value

    for (size_t idx = 0; idx < channels.size(); ++idx) {
        c->cd(idx + 1);
        const auto& data = *channels[idx].first;
        double mean = channels[idx].second.first;
        double stddev = channels[idx].second.second;
        const std::string& label = labels[idx];
        TH1D* h = makeHist(data, "h_" + label, label + " Channel;P_{#tau};Entries", mean, stddev);
        h->Draw();

        // Add blue dashed vertical line at P_tau value with label
        TLine* line = new TLine(P_tau_value, 0, P_tau_value, h->GetMaximum());
        line->SetLineColor(kBlue);
        line->SetLineStyle(2); // Dashed line
        line->SetLineWidth(2);
        line->Draw();

        // Add text label next to the line
        TLatex* txt = new TLatex(P_tau_value, h->GetMaximum() * 0.75, "P_{#tau}^{SM}");
        txt->SetTextColor(kBlue);
        txt->SetTextSize(0.05);
        txt->Draw();
    }

    const std::string output_plot = output_dir + "/segmented_summary.pdf";
    c->SaveAs(output_plot.c_str());
    std::cout << "Plot saved to: " << output_plot << std::endl;

    std::cout << ">>> All fits completed." << std::endl;
}
