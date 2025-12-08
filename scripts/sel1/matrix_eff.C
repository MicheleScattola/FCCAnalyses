#include "TString.h"
#include "TText.h"

void matrix_eff(){
    gROOT->Reset();
    
    const char* filename = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/treemaker/sel1/p8_ee_Ztautau_ecm91.root";
    const char* treename = "events";
    const int   nCat     = 6;      // categories 0,1,2,3,4,5
    const char* outdir   = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/";

    TFile *f = TFile::Open(filename, "READ");

    TTree *t = (TTree*)f->Get(treename);

    // Branches
    std::vector<int> *sel_MC_event = nullptr;
    std::vector<int> *sel_reco_event = nullptr;

    t->SetBranchAddress("sel_MC_event", &sel_MC_event);
    t->SetBranchAddress("sel_reco_event", &sel_reco_event);

    // TH2 (MC vs RECO)
    TH2D *hConf = new TH2D("hConf",
                           "Confusion Efficiency;MC;RECO",
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

    // ============================
    // normalize column
    // ============================
    for (int i = 1; i <= nCat; ++i) { 
        double col_sum = 0.0;
        for (int j = 1; j <= nCat; ++j)
            col_sum += hConf->GetBinContent(i, j);

        if (col_sum == 0.0) continue;

        for (int j = 1; j <= nCat; ++j) {
            double v = hConf->GetBinContent(i, j);
            hConf->SetBinContent(i, j, v / col_sum);   // ora 0–1
        }
    }

    // Plot 
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kStarryNight);

    TCanvas *c = new TCanvas("cConf", "Confusion matrix (Eff)", 900, 700);
    c->SetLeftMargin(0.15);
    c->SetGrid();

    hConf->GetXaxis()->SetLabelSize(0.065);  
    hConf->GetYaxis()->SetLabelSize(0.065);  

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

    hConf->GetZaxis()->SetRangeUser(0.0, 1.0);
    
    hConf->Draw("COL");

    // text 
    for (int ix = 1; ix <= nCat; ++ix) {
        double x = hConf->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= nCat; ++iy) {
            double y   = hConf->GetYaxis()->GetBinCenter(iy);
            double val = hConf->GetBinContent(ix, iy);   // 0–1

            TString label = Form("%.1f %%", val * 100.0);

            TText *t = new TText(x, y, label);
            t->SetTextAlign(22);       // centered
            t->SetTextColor(kWhite);   
            t->SetTextSize(0.035);    
            t->Draw("same");
        }
    }

    TString name_pdf = Form("%sconf_eff_sel1.pdf", outdir);

    c->SaveAs(name_pdf);

    //f->Close();
}

