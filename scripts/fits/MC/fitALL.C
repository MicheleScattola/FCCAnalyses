#include "FCCAnalyses/Fitter.h" // Adjust path if your header is elsewhere
#include <string>
#include <iostream>

void fitALL() {

    // =========================================================================
    // 1. CONFIGURATION
    // =========================================================================
    const std::string infile_data = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/treemaker/MC/p8_ee_Ztautau_ecm91.root";
    const std::string infile_templates = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/treemaker/MC/templates/templatesMC.root";
    const std::string outdir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/MC/nofilters/";
    const std::string treeName = "events";

    // Ensure output directory exists (optional system call)
    system(("mkdir -p " + outdir).c_str());

    std::cout << ">>> Starting Polarization Fits..." << std::endl;

    // =========================================================================
    // 2. TEMPLATE FITS (TFractionFitter)
    // =========================================================================
    
    // --- Rho Channel ---
    Fitter::fit(infile_data, infile_templates, outdir, 
                "TrueMC_rho.pdf",           // Output filename
                treeName, "rho_sgn",        // Data Column
                "h_template_rho_plus",      // Template +
                "h_template_rho_minus",     // Template -
                "Rho Channel (Template)",   // Title
                "#omega_{#rho}"             // X-axis label
    );

    // --- Pion Channel ---
    Fitter::fit(infile_data, infile_templates, outdir, 
                "TrueMC_pi.pdf", 
                treeName, "pi_sgn", 
                "h_template_pi_plus", 
                "h_template_pi_minus", 
                "Pion Channel (Template)", 
                "x_{#pi}"
    );

    // --- Muon Channel (Template) ---
    Fitter::fit(infile_data, infile_templates, outdir, 
                "TrueMC_mu.pdf", 
                treeName, "mu_sgn", 
                "h_template_mu_plus", 
                "h_template_mu_minus", 
                "Muon Channel (Template)", 
                "x_{#mu}"
    );

    // --- Electron Channel (Template) ---
    Fitter::fit(infile_data, infile_templates, outdir, 
                "TrueMC_el.pdf", 
                treeName, "el_sgn", 
                "h_template_el_plus", 
                "h_template_el_minus", 
                "Electron Channel (Template)", 
                "x_{e}"
    );

    // =========================================================================
    // 3. ANALYTIC FITS (Leptons Only)
    // =========================================================================
    
    // --- Muon Channel (Analytic) ---
    Fitter::fit(infile_data, outdir, 
                "TrueMC_mu_analytical.pdf", 
                treeName, "mu_sgn", 
                "Muon Channel (Analytic)", 
                "x_{#mu}",
                0.025, 1.0  // xmin, xmax
    );

    // --- Electron Channel (Analytic) ---
    Fitter::fit(infile_data, outdir, 
                "TrueMC_el_analytical.pdf", 
                treeName, "el_sgn", 
                "Electron Channel (Analytic)", 
                "x_{e}",
                0.25, 1.0  // xmin, xmax
    );

    std::cout << ">>> All fits completed." << std::endl;
}
