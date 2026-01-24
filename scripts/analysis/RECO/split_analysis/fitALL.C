#include "FCCAnalyses/Fitter.h" // Adjust path if your header is elsewhere
#include <string>
#include <iostream>

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
    Fitter::fit(infile_data, infile_templates, outdir, 
                "reco_rho.pdf",           // Output filename
                treeName, "rho_sgn",        // Data Column
                "h_template_rho_plus",      // Template +
                "h_template_rho_minus",     // Template -
                "Rho Channel (RECO)",   // Title
                "#omega_{#rho}"             // X-axis label
    );

    // --- Pion Channel ---
    Fitter::fit(infile_data, infile_templates, outdir, 
                "reco_pi.pdf", 
                treeName, "pi_sgn", 
                "h_template_pi_plus", 
                "h_template_pi_minus", 
                "Pion Channel (RECO)", 
                "x_{#pi}"
    );

    // --- Muon Channel (Template) ---
    Fitter::fit(infile_data, infile_templates, outdir, 
                "reco_mu.pdf", 
                treeName, "mu_sgn", 
                "h_template_mu_plus", 
                "h_template_mu_minus", 
                "Muon Channel (RECO)", 
                "x_{#mu}"
    );

    // --- Electron Channel (Template) ---
    Fitter::fit(infile_data, infile_templates, outdir, 
                "reco_el.pdf", 
                treeName, "el_sgn", 
                "h_template_el_plus", 
                "h_template_el_minus", 
                "Electron Channel (RECO)", 
                "x_{e}"
    );

    // =========================================================================
    // 3. ANALYTIC FITS (Leptons Only)
    // =========================================================================
    
    // --- Muon Channel (Analytic) ---
    Fitter::fit(infile_data, outdir, 
                "reco_mu_analytical.pdf", 
                treeName, "mu_sgn", 
                "Muon Channel (RECO)", 
                "x_{#mu}",
                0.05, 1.0  // xmin, xmax
    );

    // --- Electron Channel (Analytic) ---
    Fitter::fit(infile_data, outdir, 
                "reco_el_analytical.pdf", 
                treeName, "el_sgn", 
                "Electron Channel (RECO)", 
                "x_{e}",
                0.05, 1.0  // xmin, xmax
    );

    std::cout << ">>> All fits completed." << std::endl;
}
