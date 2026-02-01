#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include <sys/stat.h>

void plot() {
    const std::string input_root = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/segmented_analysis/average.root";
    const std::string output_dir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/RECO";

    TFile* fin = TFile::Open(input_root.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "[ERROR] Cannot open " << input_root << std::endl;
        return;
    }

    TH1D* h_el = (TH1D*)fin->Get("h_Electron");
    TH1D* h_mu = (TH1D*)fin->Get("h_Muon");
    TH1D* h_pi = (TH1D*)fin->Get("h_Pion");
    TH1D* h_rho = (TH1D*)fin->Get("h_Rho");

    if (!h_el || !h_mu || !h_pi || !h_rho) {
        std::cerr << "[ERROR] Missing histograms in average.root" << std::endl;
        fin->Close();
        return;
    }

    gStyle->SetOptStat(1111);
    gStyle->SetStatW(0.35);
    gStyle->SetStatH(0.25);

    TCanvas* c = new TCanvas("c", "Polarization distributions", 1200, 900);
    c->Divide(2, 2);

    const double P_tau_value = -0.14719;

    TH1D* hists[4] = {h_el, h_mu, h_pi, h_rho};
    const char* labels[4] = {"Electron", "Muon", "Pion", "Rho"};

    for (int i = 0; i < 4; ++i) {
        c->cd(i + 1);
        hists[i]->SetTitle((std::string(labels[i]) + " Channel;P_{#tau};Entries").c_str());
        hists[i]->SetLineWidth(2);
        hists[i]->SetLineColor(kBlack);
        hists[i]->Draw();

        // Fit with Gaussian using histogram mean and stddev as initial parameters
        double mean = hists[i]->GetMean();
        double stddev = hists[i]->GetStdDev();
        
        TF1* gaus = new TF1(Form("gaus_%d", i), "gaus");
        gaus->SetParameters(hists[i]->GetMaximum(), mean, stddev);
        gaus->SetLineColor(kRed);
        gaus->SetLineWidth(1);
        hists[i]->Fit(gaus, "");
        gaus->Draw("same");

        //TLine* line = new TLine(P_tau_value, 0, P_tau_value, hists[i]->GetMaximum());
        TLine* line = new TLine(P_tau_value, 0, P_tau_value, 8);
        line->SetLineColor(kBlue);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw();

        //double stddev = hists[i]->GetStdDev();
        TLatex* txt = new TLatex(P_tau_value + 0.5 * stddev, hists[i]->GetMaximum() * 0.5, "P_{#tau}^{SM}");
        txt->SetTextColor(kBlue);
        txt->SetTextSize(0.05);
        txt->Draw();
    }

    

    mkdir(output_dir.c_str(), 0755);
    const std::string output_plot = output_dir + "/segmented_summary.pdf";
    const std::string output_plot_png = output_dir + "/segmented_summary.png";
    c->SaveAs(output_plot.c_str());
    c->SaveAs(output_plot_png.c_str());
    std::cout << "Plot saved to: " << output_plot << std::endl;

    fin->Close();
}
