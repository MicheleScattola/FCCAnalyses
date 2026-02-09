#include "FCCAnalyses/Fitter.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TTree.h"

namespace {
struct WeakMixingResult {
    double sin2theta;
    double sin2theta_err;
    double gv_over_ga;
    double gv_over_ga_err;
    double A_tau;
    bool ok;
};

struct FitRecord {
    std::string label;
    double P_tau;
    double P_err;
    bool success;
    int order;
};

struct ChannelPoint {
    std::string label;
    Fitter::myFit fit;
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

std::vector<FitRecord> read_fit_results(const std::string &root_path) {
    std::vector<FitRecord> records;
    TFile in(root_path.c_str(), "READ");
    if (in.IsZombie()) {
        std::cerr << "[FitResults] Cannot open file: " << root_path
                  << std::endl;
        return records;
    }

    TTree *tree = nullptr;
    in.GetObject("fit_results", tree);
    if (!tree) {
        std::cerr << "[FitResults] Missing TTree 'fit_results' in: "
                  << root_path << std::endl;
        return records;
    }

    std::string *label = nullptr;
    double P_tau = 0.0;
    double P_err = 0.0;
    Bool_t success = false;
    int order = 0;

    tree->SetBranchAddress("label", &label);
    tree->SetBranchAddress("P_tau", &P_tau);
    tree->SetBranchAddress("P_err", &P_err);
    tree->SetBranchAddress("success", &success);
    tree->SetBranchAddress("order", &order);

    const Long64_t nentries = tree->GetEntries();
    records.reserve(static_cast<size_t>(nentries));
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        records.push_back({label ? *label : std::string(), P_tau, P_err,
                           static_cast<bool>(success), order});
    }

    std::sort(records.begin(), records.end(),
              [](const FitRecord &a, const FitRecord &b) {
                  return a.order < b.order;
              });

    return records;
}

void draw_summary_plot(const std::string &outdir,
                       const std::vector<ChannelPoint> &channels) {
    const double sm_sin2theta = 0.2315;
    const double r_sm = 1.0 - 4.0 * sm_sin2theta;
    const double A_tau_sm = 2.0 * r_sm / (1.0 + r_sm * r_sm);
    const double P_tau_sm = -A_tau_sm;

    std::vector<double> x;
    std::vector<double> xerr;
    std::vector<double> sin2;
    std::vector<double> sin2_err;
    std::vector<std::string> labels;

    for (const auto &ch : channels) {
        if (!ch.fit.success) {
            continue;
        }
        const double A_tau = -ch.fit.P_tau;
        const auto res = compute_sin2theta(A_tau, ch.fit.P_err);
        if (!res.ok) {
            continue;
        }
        const double idx = static_cast<double>(labels.size()) + 1.0;
        labels.push_back(ch.label);
        x.push_back(idx);
        xerr.push_back(0.0);
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
    const auto sin2_range = compute_range(sin2, sin2_err);
    const double axis_title_size = 0.06;
    const double axis_label_size = 0.07;

    auto draw_vertical_dividers = [](int n_bins, double ymin, double ymax) {
        std::vector<TLine *> lines;
        lines.reserve(static_cast<size_t>(n_bins));
        for (int i = 1; i < n_bins; ++i) {
            const double x = 0.5 + static_cast<double>(i);
            TLine *line = new TLine(x, ymin, x, ymax);
            line->SetLineColor(kGray + 1);
            line->SetLineStyle(3);
            line->SetLineWidth(1);
            line->Draw("SAME");
            lines.push_back(line);
        }
        return lines;
    };

    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c_summary", "Summary", 900, 500);
    gPad->SetTopMargin(0.06);
    //gPad->SetBottomMargin(0.22);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.04);

    TH1D *frame_sin2 =
        new TH1D("frame_sin2", ";;sin^{2}#theta_{W}^{eff}", n, 0.5, n + 0.5);
    for (int i = 0; i < n; ++i) {
        frame_sin2->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    }
    frame_sin2->GetYaxis()->SetRangeUser(sin2_range.first, sin2_range.second);
    frame_sin2->GetXaxis()->SetLabelSize(axis_label_size);
    frame_sin2->GetXaxis()->SetTitleSize(axis_title_size);
    frame_sin2->GetYaxis()->SetTitleSize(axis_title_size);
    frame_sin2->GetYaxis()->SetLabelSize(0.05);
    frame_sin2->GetYaxis()->SetTitleOffset(1.1);
    frame_sin2->Draw("AXIS");

    auto dividers_sin2 =
        draw_vertical_dividers(n, sin2_range.first, sin2_range.second);

    TGraphErrors *gr_sin2 = new TGraphErrors(n, x.data(), sin2.data(),
                                             xerr.data(), sin2_err.data());
    gr_sin2->SetMarkerStyle(20);
    gr_sin2->SetMarkerSize(1.0);
    gr_sin2->SetLineWidth(2);
    gr_sin2->SetLineColor(kBlack);
    gr_sin2->SetMarkerColor(kBlack);
    gr_sin2->Draw("P E1 SAME");

    TLine *line_sin2 = new TLine(0.5, sm_sin2theta, n + 0.5, sm_sin2theta);
    line_sin2->SetLineColor(kBlue);
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
    for (auto *line : dividers_sin2) {
        delete line;
    }
    delete c;
}
} // namespace

void plot_summary() {
    const std::string outdir =
        "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/RECO/split_analysis/";
    const std::string root_path = outdir + "fit_results.root";

    const auto records = read_fit_results(root_path);
    if (records.empty()) {
        std::cout << "[SummaryPlot] No records found." << std::endl;
        return;
    }

    std::vector<ChannelPoint> channels;
    channels.reserve(records.size());
    for (const auto &rec : records) {
        if (rec.label.find("Analytic") != std::string::npos) {
            continue;
        }
        Fitter::myFit fit{rec.P_tau, rec.P_err, 0.0, 0.0, rec.success,
                          rec.label};
        channels.push_back({rec.label, fit});
        print_weak_mixing(rec.label, outdir, fit);
    }

    draw_summary_plot(outdir, channels);
}
