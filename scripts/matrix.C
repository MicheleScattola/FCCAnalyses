#include "TString.h"

voidm matrix(){
	gROOT->Reset();
	
    const char* filename = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/treemaker/p8_ee_Ztautau_ecm91.root";
    const char* treename = "events";
    const int   nCat     = 6;      // categories 0,1,2,3,4,5
    const char* outdir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/";

    
    TFile *f = TFile::Open(filename, "READ");

    TTree *t = (TTree*)f->Get(treename);

    // Branches
    std::vector<int> *sel_MC_event = nullptr;
    std::vector<int> *sel_reco_event = nullptr;

    t->SetBranchAddress("sel_MC_event",        &sel_MC_event);
    t->SetBranchAddress("sel_reco_event", &sel_reco_event);

    // TH2 (MC vs RECO)
    TH2D *hConf = new TH2D("hConf",
                           "Confusion matrix;MC category;RECO category",
                           nCat, -0.5, nCat - 0.5,
                           nCat, -0.5, nCat - 0.5);

    Long64_t nEntries = t->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {

        t->GetEntry(i);

        if (!sel_MC_event || !sel_reco_event) continue;

        size_t n = sel_MC_event->size();
        if (sel_reco_event->size() != n) {
            std::cerr << "[ERROR] Mismatched sizes in sel_MC_event and sel_reco_event" << std::endl;
            continue;
        }

        for (size_t j = 0; j < n; ++j) {

            int mc   = sel_MC_event->at(j);
            int reco = sel_reco_event->at(j);

            if (mc   < 0 || mc   >= nCat) continue;
            if (reco < 0 || reco >= nCat) continue;

            hConf->Fill(mc, reco);
        }
    }

    // Plot
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kStarryNight);

    TCanvas *c = new TCanvas("cConf", "Identification matrix", 1200, 700);
    c->SetRightMargin(0.15);
    c->SetGrid();
    
    hConf->GetXaxis()->SetLabelSize(0.06);  
	hConf->GetYaxis()->SetLabelSize(0.06);  

	
	hConf->GetXaxis()->SetTickLength(0);
	hConf->GetYaxis()->SetTickLength(0);
	
	hConf->GetXaxis()->SetBinLabel(1, "other");
	hConf->GetXaxis()->SetBinLabel(2, "#mu");
	hConf->GetXaxis()->SetBinLabel(3, "e^{#pm}");
	hConf->GetXaxis()->SetBinLabel(4, "#pi^{#pm}");
	hConf->GetXaxis()->SetBinLabel(5, "#rho");
	hConf->GetXaxis()->SetBinLabel(6, "a_{1}");

	hConf->GetYaxis()->SetBinLabel(1, "other");
	hConf->GetYaxis()->SetBinLabel(2, "#mu");
	hConf->GetYaxis()->SetBinLabel(3, "e^{#pm}");
	hConf->GetYaxis()->SetBinLabel(4, "#pi^{#pm}");
	hConf->GetYaxis()->SetBinLabel(5, "#rho");
	hConf->GetYaxis()->SetBinLabel(6, "a_{1}");
	
	hConf->LabelsOption("h");  
	hConf->SetMarkerColor(kWhite);  
	hConf->SetMarkerSize(1.5);       

    hConf->Draw("COLZ TEXT");

    TString name_pdf = Form("%smatrix_id.pdf", outdir);

    //c->SaveAs(name_png);
    c->SaveAs(name_pdf);

    //f->Close();
}

