#include "FCCAnalyses/Fitter.h" // Adjust path if your header is elsewhere
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"

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

struct ChannelPoint {
    std::string label;
    const Fitter::myFit *fit;
};

void draw_summary_plot(const std::string &outdir,
                       const std::vector<ChannelPoint> &channels) {
    const double sm_sin2theta = 0.2315;
    const double r_sm = 1.0 - 4.0 * sm_sin2theta;
    const double A_tau_sm = 2.0 * r_sm / (1.0 + r_sm * r_sm);
    const double P_tau_sm = -A_tau_sm;

    std::vector<double> x;
    std::vector<double> xerr;
    std::vector<double> ptau;
    std::vector<double> ptau_err;
    std::vector<double> sin2;
    std::vector<double> sin2_err;
    std::vector<std::string> labels;

    for (const auto &ch : channels) {
        if (!ch.fit || !ch.fit->success) {
            continue;
        }
        const double A_tau = -ch.fit->P_tau;
        const auto res = compute_sin2theta(A_tau, ch.fit->P_err);
        if (!res.ok) {
            continue;
        }
        const double idx = static_cast<double>(labels.size()) + 1.0;
        labels.push_back(ch.label);
        x.push_back(idx);
        xerr.push_back(0.0);
        ptau.push_back(ch.fit->P_tau);
        ptau_err.push_back(ch.fit->P_err);
        sin2.push_back(res.sin2theta);
        sin2_err.push_back(res.sin2theta_err);
    }

    if (labels.empty()) {
        std::cout << "[SummaryPlot] No valid channels to plot." << std::endl;
        return;
    }

    auto compute_range = [](const std::vector<double> &vals,
                            const std::vector<double> &errs) {
        double ymin = vals[0] - errs[0];
        double ymax = vals[0] + errs[0];
        for (size_t i = 1; i < vals.size(); ++i) {
            ymin = std::min(ymin, vals[i] - errs[i]);
            ymax = std::max(ymax, vals[i] + errs[i]);
        }
        const double span = ymax - ymin;
        const double pad = span > 0.0 ? 0.2 * span : 0.05;
        return std::make_pair(ymin - pad, ymax + pad);
    };

    const int n = static_cast<int>(labels.size());
    const auto ptau_range = compute_range(ptau, ptau_err);
    const auto sin2_range = compute_range(sin2, sin2_err);

    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c_summary", "Summary", 900, 800);
    c->Divide(1, 2, 0.0, 0.0);

    // --- P_tau pad ---
    c->cd(1);
    gPad->SetBottomMargin(0.08);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.04);

    TH1D *frame_ptau = new TH1D("frame_ptau", ";;P_{#tau}", n, 0.5, n + 0.5);
    for (int i = 0; i < n; ++i) {
        frame_ptau->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    }
    frame_ptau->GetYaxis()->SetRangeUser(ptau_range.first, ptau_range.second);
    frame_ptau->GetXaxis()->SetLabelSize(0.04);
    frame_ptau->GetYaxis()->SetTitleSize(0.05);
    frame_ptau->GetYaxis()->SetTitleOffset(1.1);
    frame_ptau->Draw("AXIS");

    TGraphErrors *gr_ptau = new TGraphErrors(n, x.data(), ptau.data(),
                                             xerr.data(), ptau_err.data());
    gr_ptau->SetMarkerStyle(20);
    gr_ptau->SetMarkerSize(1.0);
    gr_ptau->SetLineWidth(2);
    gr_ptau->SetLineColor(kBlack);
    gr_ptau->SetMarkerColor(kBlack);
    gr_ptau->Draw("P E1 SAME");

    TLine *line_ptau = new TLine(0.5, P_tau_sm, n + 0.5, P_tau_sm);
    line_ptau->SetLineColor(kRed + 1);
    line_ptau->SetLineWidth(2);
    line_ptau->SetLineStyle(2);
    line_ptau->Draw("SAME");

    TLegend *leg_ptau = new TLegend(0.70, 0.78, 0.95, 0.93);
    leg_ptau->SetBorderSize(0);
    leg_ptau->SetFillStyle(0);
    leg_ptau->AddEntry(gr_ptau, "Fit", "ep");
    leg_ptau->AddEntry(line_ptau, "SM", "l");
    leg_ptau->Draw();

    // --- sin^2theta pad ---
    c->cd(2);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.18);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.04);

    TH1D *frame_sin2 =
        new TH1D("frame_sin2", ";;sin^{2}#theta_{W}^{eff}", n, 0.5, n + 0.5);
    for (int i = 0; i < n; ++i) {
        frame_sin2->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    }
    frame_sin2->GetYaxis()->SetRangeUser(sin2_range.first, sin2_range.second);
    frame_sin2->GetXaxis()->SetLabelSize(0.05);
    frame_sin2->GetYaxis()->SetTitleSize(0.05);
    frame_sin2->GetYaxis()->SetTitleOffset(1.1);
    frame_sin2->Draw("AXIS");

    TGraphErrors *gr_sin2 = new TGraphErrors(n, x.data(), sin2.data(),
                                             xerr.data(), sin2_err.data());
    gr_sin2->SetMarkerStyle(20);
    gr_sin2->SetMarkerSize(1.0);
    gr_sin2->SetLineWidth(2);
    gr_sin2->SetLineColor(kBlack);
    gr_sin2->SetMarkerColor(kBlack);
    gr_sin2->Draw("P E1 SAME");

    TLine *line_sin2 = new TLine(0.5, sm_sin2theta, n + 0.5, sm_sin2theta);
    line_sin2->SetLineColor(kRed + 1);
    line_sin2->SetLineWidth(2);
    line_sin2->SetLineStyle(2);
    line_sin2->Draw("SAME");

    TLegend *leg_sin2 = new TLegend(0.70, 0.78, 0.95, 0.93);
    leg_sin2->SetBorderSize(0);
    leg_sin2->SetFillStyle(0);
    leg_sin2->AddEntry(gr_sin2, "Fit", "ep");
    leg_sin2->AddEntry(line_sin2, "SM", "l");
    leg_sin2->Draw();

    const std::string pdf_path = outdir + "summary_ptau_sin2theta.pdf";
    const std::string png_path = outdir + "summary_ptau_sin2theta.png";
    c->SaveAs(pdf_path.c_str());
    c->SaveAs(png_path.c_str());

    delete leg_sin2;
    delete line_sin2;
    delete gr_sin2;
    delete frame_sin2;
    delete leg_ptau;
    delete line_ptau;
    delete gr_ptau;
    delete frame_ptau;
    delete c;
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
    print_weak_mixing("Rho (RECO)", outdir, rho_fit);

    // --- Pion Channel ---
    auto pi_fit = Fitter::fit(infile_data, infile_templates, outdir, 
                "reco_pi.pdf", 
                treeName, "pi_sgn", 
                "h_template_pi_plus", 
                "h_template_pi_minus", 
                "Pion Channel (RECO)", 
                "x_{#pi}"
    );
    print_weak_mixing("Pion (RECO)", outdir, pi_fit);

    // --- Muon Channel (Template) ---
    auto mu_fit = Fitter::fit(infile_data, infile_templates, outdir, 
                "reco_mu.pdf", 
                treeName, "mu_sgn", 
                "h_template_mu_plus", 
                "h_template_mu_minus", 
                "Muon Channel (RECO)", 
                "x_{#mu}"
    );
    print_weak_mixing("Muon (RECO)", outdir, mu_fit);

    // --- Electron Channel (Template) ---
    auto el_fit = Fitter::fit(infile_data, infile_templates, outdir, 
                "reco_el.pdf", 
                treeName, "el_sgn", 
                "h_template_el_plus", 
                "h_template_el_minus", 
                "Electron Channel (RECO)", 
                "x_{e}"
    );
    print_weak_mixing("Electron (RECO)", outdir, el_fit);

    // --- Combined Channel (Template fits) ---
    auto combined_fit = combine_channel_fits(
        "All Channels (RECO)", {&rho_fit, &pi_fit, &mu_fit, &el_fit});
    print_weak_mixing("All Channels (RECO)", outdir, combined_fit);

    // =========================================================================
    // 3. ANALYTIC FITS (Leptons Only)
    // =========================================================================
    
    // --- Muon Channel (Analytic) ---
    auto mu_ana_fit = Fitter::fit(infile_data, outdir, 
                "reco_mu_analytical.pdf", 
                treeName, "mu_sgn", 
                "Muon Channel (RECO)", 
                "x_{#mu}",
                0.05, 1.0  // xmin, xmax
    );
    print_weak_mixing("Muon (Analytic)", outdir, mu_ana_fit);

    // --- Electron Channel (Analytic) ---
    auto el_ana_fit = Fitter::fit(infile_data, outdir, 
                "reco_el_analytical.pdf", 
                treeName, "el_sgn", 
                "Electron Channel (RECO)", 
                "x_{e}",
                0.05, 1.0  // xmin, xmax
    );
    print_weak_mixing("Electron (Analytic)", outdir, el_ana_fit);

    draw_summary_plot(outdir, {
        {"Rho", &rho_fit},
        {"Pion", &pi_fit},
        {"Muon", &mu_fit},
        {"Electron", &el_fit},
        {"Combined", &combined_fit},
    });

    std::cout << ">>> All fits completed." << std::endl;
}
