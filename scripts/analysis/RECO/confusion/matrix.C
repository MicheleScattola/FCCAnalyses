#include "TString.h"

void matrix(){
	gROOT->Reset();
	
    const char* filename = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/p8_ee_Ztautau_ecm91.root";
    const char* treename = "events";
    const int   nCat     = 6;      // categories 0,1,2,3,4,5
    const char* outdir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/RECO/additional/";

    
    TFile *f = TFile::Open(filename, "READ");

    TTree *t = (TTree*)f->Get(treename);

    // Branches
    std::vector<int> *MC_event = nullptr;
    std::vector<int> *event_type_reco = nullptr;

    t->SetBranchAddress("MC_event",        &MC_event);
    t->SetBranchAddress("event_type_reco", &event_type_reco);

    // TH2 (MC vs RECO)
    TH2D *hConf = new TH2D("hConf",
                           "Confusion matrix;MC category;RECO category",
                           nCat, -0.5, nCat - 0.5,
                           nCat, -0.5, nCat - 0.5);
	
	std::vector<long> mc_counts(nCat, 0);
    Long64_t nEntries = t->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {

        t->GetEntry(i);

        if (!MC_event || !event_type_reco) continue;

        size_t n = MC_event->size();
        if (event_type_reco->size() != n) {
            std::cerr << "[ERROR] Mismatched sizes in MC_event and event_type_reco" << std::endl;
            continue;
        }

        for (size_t j = 0; j < n; ++j) {

            int mc   = MC_event->at(j);
            int reco = event_type_reco->at(j);

            if (mc   < 0 || mc   >= nCat) continue;
            if (reco < 0 || reco >= nCat) continue;

            hConf->Fill(mc, reco);
            mc_counts[mc]++;
        }
    }
    const char* labels[] = {"Other", "Mu", "El", "Pi", "Rho", "A1"};
    int total = 0;
    std::cout << "MC type counts: " << endl;
    for(int i=0; i<nCat; ++i) {
        std::cout << labels[i] << " = " << mc_counts[i] << std::endl;
        total += mc_counts[i];
    }
    std::cout << "total events = " << total/2 << endl;
    

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

