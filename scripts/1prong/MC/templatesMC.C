#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TFile.h"
#include <iostream>

void templates() {
	    
    // input
    const char* infile = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/treemaker/MC/bkg/p8_ee_Ztautau_ecm91.root";
    const char* outdir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/MC/";
    const char* outdir2 = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/treemaker/MC/bkg/";
    
    
    std::string treeName = "events"; 
    
    int nBins = 40;
    double xMin = 0.0;
    double xMax = 1.0;

    
    // Histos
    ROOT::RDataFrame df(treeName, infile);

    auto h_plus_ptr  = df.Histo1D({"h_plus",  "Template Helicity +1;x_{#pi};Events", nBins, xMin, xMax}, "pi_sgn", "w_plus");
    auto h_minus_ptr = df.Histo1D({"h_minus", "Template Helicity -1;x_{#pi};Events", nBins, xMin, xMax}, "pi_sgn", "w_minus");

    // extract histos from pointers
    TH1D *h_plus  = (TH1D*)h_plus_ptr->Clone("h_template_plus");
    TH1D *h_minus = (TH1D*)h_minus_ptr->Clone("h_template_minus");
    
    // Normalization
    h_plus->Scale(1.0 / h_plus->Integral());
    h_minus->Scale(1.0 / h_minus->Integral());

    // Plots
    gStyle->SetOptStat(0); 

    TCanvas *c = new TCanvas("c_templates", "Polarization Templates", 900, 600);
    c->SetSupportGL(true);
    c->cd();
	
    
    h_plus->SetLineColor(kBlue);
    h_plus->SetLineWidth(2);
    h_plus->SetFillColorAlpha(kBlue, 0.5); 
    //h_plus->SetFillStyle(3006);

    h_minus->SetLineColor(kRed);
    h_minus->SetLineWidth(2);
    //h_minus->SetLineStyle(2); 
    h_minus->SetFillColorAlpha(kRed, 0.5); 
    //h_minus->SetFillStyle(3007);
    
    h_plus->SetTitle("TrueMC Polarization Templates;x (re-weighted);Probability Density");
    

    
    double max_y = std::max(h_plus->GetMaximum(), h_minus->GetMaximum());
    h_plus->SetMaximum(max_y * 1.2); // +20% margin
    h_plus->SetMinimum(0.);

    h_plus->Draw("HIST");
    h_minus->Draw("HIST SAME");

    // Legend
    TLegend *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(h_plus,  "Helicity +1", "l");
    leg->AddEntry(h_minus, "Helicity -1", "l");
    leg->Draw();

    // save canvas
    //c->SaveAs(TString(outdir) + "templates.png");
    c->SaveAs(TString(outdir) + "templates.pdf");

    // SAVING HISTOS FOR FINAL FIT
    TString rootOutName = TString(outdir2) + "templatesMC.root";
    TFile *fOut = new TFile(rootOutName, "RECREATE");
    h_plus->Write();
    h_minus->Write();
    fOut->Close();

    std::cout << "\n[INFO] Templates created and saved in: " << rootOutName << std::endl;
    std::cout << "       - h_template_plus" << std::endl;
    std::cout << "       - h_template_minus" << std::endl;
}
