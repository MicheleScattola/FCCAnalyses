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

double skewed_gaus(double *x, double *p)
{
  double mean     = p[0];   // mean
  double sigma    = p[1];   // sigma
  double skewness  = p[2];   // skewness
  double amplitude = p[3];   // amplitude

  double argument = (x[0] - mean)/ sigma;

  double Norm_gaus = TMath::Gaus(x[0], mean, sigma, true);
  // cumulative distribution function
  double CDF       = 0.5 * (1 + erf(skewness* argument/sqrt(2)) );

  return amplitude * Norm_gaus * CDF;
}

void Fit_skewed(TH1 *hist) {
        
    if (!hist) {
        std::cerr << "[ERROR] Fit_skewed: hist is null\n";
        return;
    }

    
    TH1 *h = (TH1*)hist->Clone(Form("%s_clone", hist->GetName()));
    
    // normalization
    double integral = h->Integral("width"); 
    int N = h->Integral();
    
    if (integral <= 0) {
        std::cerr << "[ERROR] Fit_skewed: hist integral <= 0\n";
        delete h;
        return;
    } else {
    	std::cout << "\n =========================\n" << endl;
        std::cout << "[INFO] Fitting " << hist->GetName() << " with " << N << " entries." << std::endl;
    }
    
    h->Scale(1.0 / integral);

    
    double xmin = 0.05;
    double xmax = 0.95;
    
    // Polynomial FIT
    TF1 *f = new TF1("bkg_shape", skewed_gaus, xmin, xmax,4);
    
    f->SetLineColor(kRed);
    f->SetLineWidth(2);
    
    // initialize values
	f->SetParNames("x_{0}","#sigma","#alpha","A");
	f->SetParameter(0, 0.2); 
    f->SetParameter(1, 0.5);  
    f->SetParameter(2, 12);   
    f->SetParameter(3, 4); 

    std::cout << "Fitting skewed gaussian between " << xmin << " and " << xmax << std::endl;      

    gStyle->SetOptStat(0);      
    gStyle->SetOptFit(1111);    
    
    // plot
    TCanvas *c = new TCanvas(Form("canvas_%s", hist->GetName()), "bkg skewed gaus fit", 800,600);
    c->cd();
    
    h->SetTitle("Background fit (skewed gaussian);x;Normalized Events");
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.7);
    h->SetLineColor(kBlack);
    h->SetLineWidth(1);
    
    h->SetStats(1);
	h->Draw("HIST E");
	h->Fit(f, "R SAME");
	f->Draw("SAME");
    
    std::cout << "\n=========================\n" << std::endl;
    
    c->Modified();
    c->Update();

    const char* outdir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/plots/optimal/";
    
    
    c->SaveAs(TString(outdir) + "piBKG_skew.png");
    c->SaveAs(TString(outdir) + "piBKG_skew.pdf");
    
    // delete h; 
    // delete c;
    // delete f;
}


void bkg_skewed() {

    const char* infile = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/histmaker/optimal/p8_ee_Ztautau_ecm91.root";

    TFile* f = TFile::Open(infile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open file: " << infile << std::endl;
        return;
    }

   
    TH1F* h = (TH1F*)f->Get("pi_bkg");

    if (!h) {
        std::cerr << "[WARN] Histogram pi_bkg not found in file" << std::endl;
        f->Close(); // Chiudi il file se fallisce
        return;
    }

    Fit_skewed(h);

    
    // f->Close(); 
}
