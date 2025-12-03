#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TMath.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "ROOT/RDataFrame.hxx"
#include <iostream>

// Global pointers for chi2 function
static TH1D *g_data  = nullptr;
static TH1D *g_plus  = nullptr;
static TH1D *g_minus = nullptr;

double chi2_function(const double *par) {
    //
    // par[0] = f_plus
    //
    double f = par[0];
    if (f < 0 || f > 1) return 1e30;  // enforce physical boundary

    double chi2 = 0.0;

    int nb = g_data->GetNbinsX();
    for (int i = 1; i <= nb; i++) {

        double d = g_data->GetBinContent(i);
        double e = g_data->GetBinError(i);

        // skip bins with zero data uncertainty
        if (e <= 0) continue;

        // model = linear combination of templates
        double m = f * g_plus->GetBinContent(i) +
                   (1 - f) * g_minus->GetBinContent(i);

        chi2 += (d - m) * (d - m) / (e * e);
    }

    return chi2;
}



void fitFraction_chi2() {

    const char* infile_data = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/treemaker/optimal/p8_ee_Ztautau_ecm91.root";
    const char* infile_templates = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/treemaker/optimal/bkg/templates.root";
    std::string treeName = "events";

    // --- LOAD TEMPLATES ---
    TFile *fTemp = TFile::Open(infile_templates);
    TH1D *h_plus  = (TH1D*)fTemp->Get("h_template_plus");
    TH1D *h_minus = (TH1D*)fTemp->Get("h_template_minus");
    h_plus->SetDirectory(0);
    h_minus->SetDirectory(0);
    fTemp->Close();

    // --- LOAD DATA ---
    ROOT::RDataFrame df(treeName, infile_data);
    auto h_data_ptr = df.Histo1D(
        {"h_data", "Data", h_plus->GetNbinsX(),
         h_plus->GetXaxis()->GetXmin(),
         h_plus->GetXaxis()->GetXmax()},
        "pi_sgn"
    );
    TH1D *h_data = (TH1D*)h_data_ptr->Clone("h_data");
    h_data->Sumw2();

    // --- Set globals for chi2 ---
    g_data  = h_data;
    g_plus  = h_plus;
    g_minus = h_minus;

    // --- MINIMIZER ---
    ROOT::Minuit2::Minuit2Minimizer min(ROOT::Minuit2::kMigrad);

    ROOT::Math::Functor f(&chi2_function, 1);
    min.SetFunction(f);

    // parameter: f_plus
    min.SetVariable(0, "f_plus", 0.5, 0.01);
    min.SetVariableLimits(0, 0.0, 1.0);

    min.Minimize();

    double f_plus  = min.X()[0];
    double err_f   = min.Errors()[0];
    double f_minus = 1.0 - f_plus;

    // Polarization
    double P     = f_plus - f_minus;
    double err_P = 2 * err_f;

    std::cout << "=== Fit results (chi2) ===" << std::endl;
    std::cout << "f_plus  = " << f_plus  << " +/- " << err_f   << std::endl;
    std::cout << "f_minus = " << f_minus << " +/- " << err_f   << std::endl;
    std::cout << "P       = " << P       << " +/- " << err_P   << std::endl;
    std::cout << "Chi2    = " << min.MinValue() << std::endl;

    // --- Produce fitted model histogram ---
    TH1D *h_fit = (TH1D*)h_plus->Clone("h_fit");
    for (int i = 1; i <= h_fit->GetNbinsX(); i++) {
        double model = f_plus * h_plus->GetBinContent(i)
                     + (1 - f_plus) * h_minus->GetBinContent(i);
        h_fit->SetBinContent(i, model);
        h_fit->SetBinError(i, 0);
    }
    h_fit->SetLineColor(kRed);

    // --- PLOT ---
    TCanvas *c = new TCanvas("c","c",800,600);
    h_data->SetMarkerStyle(20);
    h_data->Draw("E");
    h_fit->Draw("HIST SAME");
}

