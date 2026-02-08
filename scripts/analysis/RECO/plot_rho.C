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

void plot_rho() {
	const std::string infile =
			"/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/treemaker/p8_ee_Ztautau_ecm91.root";
	const std::string outdir = get_dirname(infile);

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	//gStyle->SetTitleSize(0.06, "t");

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

	if (!tree->GetBranch("rho_true") || !tree->GetBranch("rho_fake")) {
		std::cerr << "[plot_rho] Missing branches rho_true and/or rho_fake." << std::endl;
		file->Close();
		return;
	}

	TH1D *h_true = new TH1D("h_rho_true", "True #rho VS a_{1};#omega_{#rho};Events", 40, -1.0, 1.0);
	TH1D *h_fake = new TH1D("h_rho_fake", "#rho_{fake};#omega_{#rho};Events", 40, -1.0, 1.0);

	h_true->GetXaxis()->SetTitleSize(0.055);
	h_true->GetYaxis()->SetTitleSize(0.045);
	h_fake->GetXaxis()->SetTitleSize(0.055);
	h_fake->GetYaxis()->SetTitleSize(0.045);

	tree->Draw("rho_true>>h_rho_true", "", "goff");
	tree->Draw("rho_fake>>h_rho_fake", "", "goff");

	const double true_area = h_true->Integral();
	const double fake_area = h_fake->Integral();
	if (true_area > 0.0 && fake_area > 0.0) {
		h_fake->Scale((true_area / 3.0) / fake_area);
	}

	TCanvas *c = new TCanvas("c_rho", "Rho", 800, 700);
	gPad->SetTopMargin(0.08);
	gPad->SetBottomMargin(0.14);

	h_true->SetLineColor(kBlue + 1);
	h_true->SetLineWidth(2);
	h_true->Draw("HIST");

	h_fake->SetLineColor(kRed + 1);
	h_fake->SetLineWidth(2);
	h_fake->Draw("HIST SAME");

	TLegend *leg = new TLegend(0.55, 0.65, 0.88, 0.88);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->AddEntry(h_true, "#rho_{true}", "l");
	leg->AddEntry(h_fake, "a_{1} BKG (re-scaled)", "l");
	leg->Draw();

	const std::string pdf_out = outdir + "rho_true_fake.pdf";
	const std::string png_out = outdir + "rho_true_fake.png";
	c->SaveAs(pdf_out.c_str());
	c->SaveAs(png_out.c_str());

	delete leg;
	delete h_fake;
	delete h_true;
	delete c;
	file->Close();
}
