#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TMath.h"

void thrust() {
   
    std::string inputDir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/treemaker/";
    std::string fileName = "p8_ee_Ztautau_ecm91"; 
    std::string inputFile = inputDir + fileName + ".root";
    
    std::string outputDir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/";
    std::string outputFile = outputDir + "thrust.pdf";
    std::string treeName = "events"; 

    gStyle->SetOptStat(0); 
    gStyle->SetOptTitle(0); 

    // load histo with dataframe
    ROOT::RDataFrame df(treeName, inputFile);
    // CosTheta [-1, 1]
    auto h_costheta_ptr = df.Histo1D({"h_costheta", ";Thrust cos(#theta);Events", 100, -1.0, 1.0}, "RP_thrustcostheta");
    // Phi [0, pi]
    auto h_phi_ptr = df.Histo1D({"h_phi", ";Thrust #phi;Events", 100, 0.0, TMath::Pi()}, "RP_thrustphi");

    TH1D *h_costheta = (TH1D*)h_costheta_ptr->Clone();
    TH1D *h_phi = (TH1D*)h_phi_ptr->Clone();


    h_costheta->SetLineColor(kBlue + 1);
    h_costheta->SetLineWidth(2);
    h_costheta->SetFillColorAlpha(kBlue + 1, 0.3);

    h_phi->SetLineColor(kRed + 1);
    h_phi->SetLineWidth(2);
    h_phi->SetFillColorAlpha(kRed + 1, 0.3);

    // axis
    h_costheta->SetMinimum(0); 
    //h_costheta->SetMaximum(h_costheta->GetMaximum() * 1.35); 

    h_phi->SetMinimum(0); 
    h_phi->SetMaximum(h_phi->GetMaximum() * 1.20); 

 
    float labelSize = 0.045;
    float titleSize = 0.05;
    
    h_costheta->GetYaxis()->SetLabelSize(labelSize);
    h_costheta->GetYaxis()->SetTitleSize(titleSize);
    h_costheta->GetXaxis()->SetLabelSize(labelSize);
    h_costheta->GetXaxis()->SetTitleSize(titleSize);

    h_phi->GetYaxis()->SetLabelSize(labelSize);
    h_phi->GetYaxis()->SetTitleSize(titleSize);
    h_phi->GetXaxis()->SetLabelSize(labelSize);
    h_phi->GetXaxis()->SetTitleSize(titleSize);

    // canvas
    TCanvas *c = new TCanvas("c", "Thrust Variables", 900, 800);
    // Divisione standard 1x2 (con margini di default tra i pad)
    c->Divide(1, 2); 

    // COSTHETA
    c->cd(1);
    
    h_costheta->Draw("HIST");

    // center legend
    TLegend *leg1 = new TLegend(0.35, 0.75, 0.65, 0.92);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.06);
    leg1->SetTextAlign(22); 
    leg1->AddEntry(h_costheta, "Thrust cos(#theta)", "f");
    leg1->Draw();

    // PHI
    c->cd(2);
    
    h_phi->Draw("HIST");

    // Legenda centrata in alto
    TLegend *leg2 = new TLegend(0.35, 0.80, 0.65, 0.92);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.06);
    leg2->SetTextAlign(22);
    leg2->AddEntry(h_phi, "Thrust #phi", "f");
    leg2->Draw();

    // save canvas
    gSystem->Exec(("mkdir -p " + outputDir).c_str());
    c->SaveAs(outputFile.c_str());
    
}
