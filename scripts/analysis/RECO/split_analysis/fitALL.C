#include "FCCAnalyses/Fitter.h" // Adjust path if your header is elsewhere
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"

namespace {
struct FitRecord {
    std::string label;
    double P_tau;
    double P_err;
    bool success;
    int order;
};

Fitter::myFit combine_channel_fits(
    const std::string &label,
    const std::vector<const Fitter::myFit *> &fits) {
    Fitter::myFit result{0.0, 0.0, 0.0, 0.0, false, "Combined"};
    double sum_w = 0.0;
    double sum_wA = 0.0;
    int used = 0;

    for (const auto *fit : fits) {
        if (!fit || !fit->success) {
            continue;
        }
        if (!std::isfinite(fit->P_tau) || !std::isfinite(fit->P_err) ||
            fit->P_err <= 0.0) {
            continue;
        }
        const double A_tau = -fit->P_tau;
        const double w = 1.0 / (fit->P_err * fit->P_err);
        sum_w += w;
        sum_wA += w * A_tau;
        ++used;
    }

    if (sum_w <= 0.0) {
        std::cout << "[Combine] " << label << ": no valid fits."
                  << std::endl;
        return result;
    }

    const double A_tau = sum_wA / sum_w;
    const double A_err = std::sqrt(1.0 / sum_w);
    result.P_tau = -A_tau;
    result.P_err = A_err;
    result.success = true;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "[Combine] " << label << " | A_tau = " << A_tau
              << " +/- " << A_err << " (used " << used << " channels)"
              << std::endl;
    return result;
}

void write_fit_results(const std::string &outdir,
                       const std::vector<FitRecord> &records) {
    const std::string root_path = outdir + "fit_results.root";
    TFile out(root_path.c_str(), "RECREATE");
    if (out.IsZombie()) {
        std::cerr << "[FitResults] Cannot create file: " << root_path
                  << std::endl;
        return;
    }

    TTree tree("fit_results", "fit results");
    std::string label;
    double P_tau = 0.0;
    double P_err = 0.0;
    Bool_t success = false;
    int order = 0;

    tree.Branch("label", &label);
    tree.Branch("P_tau", &P_tau);
    tree.Branch("P_err", &P_err);
    tree.Branch("success", &success);
    tree.Branch("order", &order);

    for (const auto &rec : records) {
        label = rec.label;
        P_tau = rec.P_tau;
        P_err = rec.P_err;
        success = rec.success;
        order = rec.order;
        tree.Fill();
    }

    tree.Write();
    out.Close();
}
} // namespace

void fitALL() {

    // =========================================================================
    // 1. CONFIGURATION
    // =========================================================================
    const std::string infile_data = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/split_analysis/output_merged.root";
    const std::string infile_templates = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/split_analysis/templates_histograms.root";
    const std::string outdir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/RECO/split_analysis/";
    const std::string treeName = "events";

    // Ensure output directory exists (optional system call)
    system(("mkdir -p " + outdir).c_str());

    std::cout << ">>> Starting Polarization Fits..." << std::endl;

    // =========================================================================
    // 2. TEMPLATE FITS 
    // =========================================================================
    
    // --- Rho Channel ---
    auto rho_fit = Fitter::fit(infile_data, infile_templates, outdir, 
                "reco_rho.pdf",           // Output filename
                treeName, "rho_sgn",        // Data Column
                "h_template_rho_plus",      // Template +
                "h_template_rho_minus",     // Template -
                "Rho Channel (RECO)",   // Title
                "#omega_{#rho}"             // X-axis label
    );
    std::vector<FitRecord> records;
    int order = 0;
    auto push_record = [&](const std::string &label,
                           const Fitter::myFit &fit) {
        records.push_back({label, fit.P_tau, fit.P_err, fit.success, order});
        ++order;
    };
    push_record("Rho", rho_fit);

    // --- Pion Channel ---
    auto pi_fit = Fitter::fit(infile_data, infile_templates, outdir, 
                "reco_pi.pdf", 
                treeName, "pi_sgn", 
                "h_template_pi_plus", 
                "h_template_pi_minus", 
                "Pion Channel (RECO)", 
                "x_{#pi}"
    );
    push_record("Pion", pi_fit);

    // --- Muon Channel (Template) ---
    auto mu_fit = Fitter::fit(infile_data, infile_templates, outdir, 
                "reco_mu.pdf", 
                treeName, "mu_sgn", 
                "h_template_mu_plus", 
                "h_template_mu_minus", 
                "Muon Channel (RECO)", 
                "x_{#mu}"
    );
    push_record("Muon", mu_fit);

    // --- Electron Channel (Template) ---
    auto el_fit = Fitter::fit(infile_data, infile_templates, outdir, 
                "reco_el.pdf", 
                treeName, "el_sgn", 
                "h_template_el_plus", 
                "h_template_el_minus", 
                "Electron Channel (RECO)", 
                "x_{e}"
    );
    push_record("Electron", el_fit);

    // --- Combined Channel (Template fits) ---
    auto combined_fit = combine_channel_fits(
        "All Channels (RECO)", {&rho_fit, &pi_fit, &mu_fit, &el_fit});
    push_record("Combined", combined_fit);

    write_fit_results(outdir, records);

    std::cout << ">>> All fits completed." << std::endl;
}
