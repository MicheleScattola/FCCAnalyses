#include "FCCAnalyses/Fitter.h" // Adjust path if your header is elsewhere
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace {
struct WeakMixingResult {
    double sin2theta;
    double sin2theta_err;
    double gv_over_ga;
    double gv_over_ga_err;
    double A_tau;
    bool ok;
};

WeakMixingResult compute_sin2theta(double A_tau, double A_tau_err) {
    WeakMixingResult out{0.0, 0.0, 0.0, 0.0, A_tau, false};
    const double absA = std::fabs(A_tau);
    if (absA == 0.0 || absA >= 1.0) {
        return out;
    }

    const double s = std::sqrt(1.0 - A_tau * A_tau);
    const double r = (1.0 - s) / A_tau; // gv/ga, physical small-r branch
    const double sin2theta = (1.0 - r) / 4.0;

    const double dr_dA = (1.0 / s) - (r / A_tau);
    const double ds2_dA = -0.25 * dr_dA;
    const double gv_over_ga_err = std::fabs(dr_dA) * A_tau_err;
    const double sin2theta_err = std::fabs(ds2_dA) * A_tau_err;

    out.sin2theta = sin2theta;
    out.sin2theta_err = sin2theta_err;
    out.gv_over_ga = r;
    out.gv_over_ga_err = gv_over_ga_err;
    out.ok = true;
    return out;
}

void print_weak_mixing(const std::string &label, const std::string &outdir,
                                             const Fitter::myFit &fit) {
    if (!fit.success) {
        std::cout << "[WeakMixing] " << label << ": fit failed." << std::endl;
        return;
    }

    const double A_tau = -fit.P_tau;
    const auto res = compute_sin2theta(A_tau, fit.P_err);
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "[WeakMixing] " << label
              << " | A_tau = " << A_tau << " +/- " << fit.P_err;

    if (!res.ok) {
        std::cout << " | sin^2theta: invalid A_tau" << std::endl;
        return;
    }

    std::cout << " | gv/ga = " << res.gv_over_ga
              << " +/- " << res.gv_over_ga_err
              << " | sin^2theta = " << res.sin2theta
              << " +/- " << res.sin2theta_err << std::endl;

      if (outdir.empty()) {
          return;
      }

      const std::string txt_path = outdir + "weak_mixing_results.txt";
      const bool write_header = !std::ifstream(txt_path).good();
      std::ofstream out(txt_path, std::ios::app);
      if (!out.good()) {
          std::cerr << "[WeakMixing] Warning: cannot open file: " << txt_path
                    << std::endl;
          return;
      }

      if (write_header) {
          out << "label\tA_tau\tA_tau_err\tgv_over_ga\tgv_over_ga_err\tsin2theta\tsin2theta_err"
              << std::endl;
      }

    out << label << "\t" << A_tau << "\t" << fit.P_err << "\t"
          << res.gv_over_ga << "\t" << res.gv_over_ga_err << "\t"
          << res.sin2theta << "\t" << res.sin2theta_err << std::endl;
}
} // namespace

void fitALL() {

    // =========================================================================
    // 1. CONFIGURATION
    // =========================================================================
    const std::string infile_data = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/p8_ee_Ztautau_ecm91.root";
    const std::string infile_templates = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/templates/templates.root";
    const std::string outdir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/RECO/fits/";
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
                "Rho Channel (Template)",   // Title
                "#omega_{#rho}"             // X-axis label
    );
    print_weak_mixing("Rho (Template)", outdir, rho_fit);

    // --- Pion Channel ---
    auto pi_fit = Fitter::fit(infile_data, infile_templates, outdir, 
                "reco_pi.pdf", 
                treeName, "pi_sgn", 
                "h_template_pi_plus", 
                "h_template_pi_minus", 
                "Pion Channel (Template)", 
                "x_{#pi}"
    );
    print_weak_mixing("Pion (Template)", outdir, pi_fit);

    // --- Muon Channel (Template) ---
    auto mu_fit = Fitter::fit(infile_data, infile_templates, outdir, 
                "reco_mu.pdf", 
                treeName, "mu_sgn", 
                "h_template_mu_plus", 
                "h_template_mu_minus", 
                "Muon Channel (Template)", 
                "x_{#mu}"
    );
    print_weak_mixing("Muon (Template)", outdir, mu_fit);

    // --- Electron Channel (Template) ---
    auto el_fit = Fitter::fit(infile_data, infile_templates, outdir, 
                "reco_el.pdf", 
                treeName, "el_sgn", 
                "h_template_el_plus", 
                "h_template_el_minus", 
                "Electron Channel (Template)", 
                "x_{e}"
    );
    print_weak_mixing("Electron (Template)", outdir, el_fit);

    // =========================================================================
    // 3. ANALYTIC FITS (Leptons Only)
    // =========================================================================
    
    // --- Muon Channel (Analytic) ---
    auto mu_ana_fit = Fitter::fit(infile_data, outdir, 
                "reco_mu_analytical.pdf", 
                treeName, "mu_sgn", 
                "Muon Channel (Analytic)", 
                "x_{#mu}",
                0.05, 1.0  // xmin, xmax
    );
    print_weak_mixing("Muon (Analytic)", outdir, mu_ana_fit);

    // --- Electron Channel (Analytic) ---
    auto el_ana_fit = Fitter::fit(infile_data, outdir, 
                "reco_el_analytical.pdf", 
                treeName, "el_sgn", 
                "Electron Channel (Analytic)", 
                "x_{e}",
                0.05, 1.0  // xmin, xmax
    );
    print_weak_mixing("Electron (Analytic)", outdir, el_ana_fit);

    std::cout << ">>> All fits completed." << std::endl;
}
