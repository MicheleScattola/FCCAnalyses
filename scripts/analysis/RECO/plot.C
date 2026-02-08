#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TKey.h"

namespace {
TTree *get_tree(TFile *file) {
	if (!file) {
		return nullptr;
	}
	TTree *tree = dynamic_cast<TTree *>(file->Get("events"));
	if (tree) {
		return tree;
	}
	TIter next(file->GetListOfKeys());
	TKey *key = nullptr;
	while ((key = static_cast<TKey *>(next()))) {
		TObject *obj = key->ReadObj();
		tree = dynamic_cast<TTree *>(obj);
		if (tree) {
			return tree;
		}
	}
	return nullptr;
}

std::string get_dirname(const std::string &path) {
	const std::string::size_type pos = path.find_last_of('/');
	if (pos == std::string::npos) {
		return "./";
	}
	return path.substr(0, pos + 1);
}
} // namespace

void plot() {
	const std::string infile =
			"/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/treemaker/p8_ee_Ztautau_ecm91.root";
	const std::string outdir = get_dirname(infile);

	const double expected_mean = -0.011;
	const double expected_sigma = 0.05;
	const double fit_min = expected_mean - 1.0 * expected_sigma;
	const double fit_max = expected_mean + 1.0 * expected_sigma;

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1111);
	gStyle->SetTitleSize(0.06, "t");
	gStyle->SetStatX(0.50);
	gStyle->SetStatY(0.90);

	TFile *file = TFile::Open(infile.c_str(), "READ");
	if (!file || file->IsZombie()) {
		std::cerr << "[plot_sys] Cannot open file: " << infile << std::endl;
		return;
	}

	TTree *tree = get_tree(file);
	if (!tree) {
		std::cerr << "[plot_sys] No TTree found in file: " << infile << std::endl;
		file->Close();
		return;
	}

	if (!tree->GetBranch("e_sys") || !tree->GetBranch("mu_sys")) {
		std::cerr << "[plot_sys] Missing branches e_sys and/or mu_sys." << std::endl;
		file->Close();
		return;
	}

	TH1D *h_e = new TH1D("h_e_sys", "e_{sys};(#Delta x/x)_e;Events", 80, -0.015, 0.007);
	TH1D *h_mu = new TH1D("h_mu_sys", "#mu_{sys};(#Delta x/x)_#mu;Events", 80, -0.015, 0.007);

	h_e->GetXaxis()->SetTitleSize(0.055);
	h_e->GetYaxis()->SetTitleSize(0.055);
	h_mu->GetXaxis()->SetTitleSize(0.055);
	h_mu->GetYaxis()->SetTitleSize(0.055);

	tree->Draw("e_sys>>h_e_sys", "", "goff");
	tree->Draw("mu_sys>>h_mu_sys", "", "goff");

	TF1 *gaus_e = new TF1("gaus_e", "gaus", fit_min, fit_max);
	gaus_e->SetParameters(h_e->GetMaximum(), expected_mean, expected_sigma);
	h_e->Fit(gaus_e, "Q");

	TF1 *gaus_mu = new TF1("gaus_mu", "gaus", fit_min, fit_max);
	gaus_mu->SetParameters(h_mu->GetMaximum(), expected_mean, expected_sigma);
	h_mu->Fit(gaus_mu, "Q");

	TCanvas *c = new TCanvas("c_sys", "Systematics", 800, 900);
	c->Divide(1, 2, 0.0, 0.0);

	c->cd(1);
	gPad->SetTopMargin(0.08);
	gPad->SetBottomMargin(0.20);
	h_e->SetLineColor(kBlue + 1);
	h_e->SetLineWidth(2);
	h_e->Draw("HIST");
	gaus_e->SetLineColor(kRed + 1);
	gaus_e->SetLineWidth(2);
	gaus_e->Draw("SAME");

	TLegend *leg_e = new TLegend(0.60, 0.72, 0.89, 0.88);
	leg_e->SetBorderSize(0);
	leg_e->SetFillStyle(0);
	leg_e->AddEntry(h_e, "e_{sys}", "l");
	leg_e->AddEntry(gaus_e, "Gaussian fit", "l");
	//leg_e->Draw();

	c->cd(2);
	gPad->SetTopMargin(0.08);
	gPad->SetBottomMargin(0.20);
	h_mu->SetLineColor(kBlue + 1);
	h_mu->SetLineWidth(2);
	h_mu->Draw("HIST");
	gaus_mu->SetLineColor(kRed + 1);
	gaus_mu->SetLineWidth(2);
	gaus_mu->Draw("SAME");

	TLegend *leg_mu = new TLegend(0.60, 0.72, 0.89, 0.88);
	leg_mu->SetBorderSize(0);
	leg_mu->SetFillStyle(0);
	leg_mu->AddEntry(h_mu, "#mu_{sys}", "l");
	leg_mu->AddEntry(gaus_mu, "Gaussian fit", "l");
	//leg_mu->Draw();

	const std::string pdf_out = outdir + "sys_gaussian_fits.pdf";
	const std::string png_out = outdir + "sys_gaussian_fits.png";
	c->SaveAs(pdf_out.c_str());
	c->SaveAs(png_out.c_str());

	delete leg_mu;
	delete leg_e;
	delete gaus_mu;
	delete gaus_e;
	delete h_mu;
	delete h_e;
	delete c;
	file->Close();
}
