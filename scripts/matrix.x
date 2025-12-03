#include "TString.h"

{
	gROOT->Reset();
	
    const char* filename = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/treemaker/optimal/p8_ee_Ztautau_ecm91.root";
    const char* treename = "events";
    const int   nCat     = 6;      // categories 0,1,2,3,4,5
    const char* outdir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/plots/confusion";

    
    TFile *f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) {
        printf("Errore: impossibile aprire %s\n", filename);
        return;
    }

    TTree *t = (TTree*)f->Get(treename);
    if (!t) {
        printf("Error: TTree %s not found\n", treename);
        f->Close();
        return;
    }

    // Branches
    std::vector<int> *MC_event        = nullptr;
    std::vector<int> *event_type_reco = nullptr;

    t->SetBranchAddress("MC_event",        &MC_event);
    t->SetBranchAddress("event_type_reco", &event_type_reco);

    // TH2 (MC vs RECO)
    TH2D *hConf = new TH2D("hConf",
                           "Confusion matrix;MC category;RECO category",
                           nCat, -0.5, nCat - 0.5,
                           nCat, -0.5, nCat - 0.5);

    Long64_t nEntries = t->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {

        t->GetEntry(i);

        if (!MC_event || !event_type_reco) continue;

        size_t n = MC_event->size();
        if (event_type_reco->size() != n) {
            printf("Warning entry %lld: size MC=%zu RECO=%zu\n",
                   i, MC_event->size(), event_type_reco->size());
            n = std::min(MC_event->size(), event_type_reco->size());
        }

        for (size_t j = 0; j < n; ++j) {

            int mc   = MC_event->at(j);
            int reco = event_type_reco->at(j);

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

	// opzionale: togliere anche i tick
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
	hConf->SetMarkerColor(kWhite);   // colore del testo nelle celle
	hConf->SetMarkerSize(1.5);       // opzionale: testo più grande

    hConf->Draw("COLZ TEXT");

    TString name_png = Form("%s/matrix_id.png", outdir);
    TString name_pdf = Form("%s/matrix_id.pdf", outdir);

    c->SaveAs(name_png);
    c->SaveAs(name_pdf);

    //f->Close();
}

