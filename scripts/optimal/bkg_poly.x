#include <iostream>
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TString.h"

//================================================
// Funzione di Fitting
//================================================
void Fit_bkg_poly(TH1 *hist) {
        
    if (!hist) {
        cerr << "[ERROR] Fit_bkg_poly: hist is null\n";
        return;
    }

    
    TH1 *h = (TH1*)hist->Clone(Form("%s_clone", hist->GetName()));
    h->Sumw2();
    double N = h->Integral();
    
	cout << "\n =========================\n" << endl;
    cout << "[INFO] Fitting " << hist->GetName() << " with " << N << " entries." << endl;

    double xmin = 0.05;
    double xmax = 0.95;
    
    // Polynomial FIT
    TF1 *f = new TF1("bkg_shape", "pol7", xmin, xmax);
    
    f->SetLineColor(kRed);
    f->SetLineWidth(2);

    cout << "Fitting polynomial of degree 7 between " << xmin << " and " << xmax << endl;      

    gStyle->SetOptStat(0);      
    gStyle->SetOptFit(1111);    
    
    // plot
    TCanvas *c = new TCanvas(Form("canvas_%s", hist->GetName()), "bkg poly fit", 900,600);
    c->cd();
    
    h->SetTitle("Background fit (Pol7);x;Events");
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.7);
    h->SetLineColor(kBlack);
    h->SetLineWidth(1);
    
    
    // "M" = (Minuit improvement) - opzionale con polN 
    h->SetStats(1);
    h->Draw("HIST E");
    //h->Fit(f, "R  +"); 
    h->Fit(f, "R M +"); 
    f->Draw("SAME");
    
    cout << "\n=========================\n" << endl;
    
    c->Modified();
    c->Update();

    const char* outdir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/plots/optimal/";
    
    // saving shape into root file
    TFile *fOut = new TFile("bkg_parameters.root", "RECREATE");
    cout << "[INFO] saving background shape into file 'bkg_parameters.root' \n" << endl;
    f->Write();
    fOut->Close();
    // saving canvas
    c->SaveAs(TString(outdir) + "piBKG_poly.png");
    c->SaveAs(TString(outdir) + "piBKG_poly.pdf");
    
    // delete h; 
    // delete c;
    // delete f;
}


void bkg_poly() {

    const char* infile = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/histmaker/optimal/p8_ee_Ztautau_ecm91.root";

    TFile* f = TFile::Open(infile, "READ");
    if (!f || f->IsZombie()) {
        cerr << "ERROR: cannot open file: " << infile << endl;
        return;
    }

   
    TH1F* h = (TH1F*)f->Get("pi_bkg");

    if (!h) {
        cerr << "[WARN] Histogram pi_bkg not found in file" << endl;
        f->Close(); // Chiudi il file se fallisce
        return;
    }

    Fit_bkg_poly(h);

    
    f->Close(); 
}
