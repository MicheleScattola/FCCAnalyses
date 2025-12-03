#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TMath.h"

void thrust() {
    // ================= CONFIGURAZIONE =================
    std::string inputDir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/treemaker/";
    std::string fileName = "p8_ee_Ztautau_ecm91"; 
    std::string inputFile = inputDir + fileName + ".root";
    
    std::string outputDir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/plots/";
    std::string outputFile = outputDir + "thrust_final.pdf";
    std::string treeName = "events"; 
    // ==================================================

    // Rimuove il titolo superiore del plot e il box statistico
    gStyle->SetOptStat(0); 
    gStyle->SetOptTitle(0); 

    // 1. Apertura File e Dataframe
    ROOT::RDataFrame df(treeName, inputFile);

    // 2. Definizione istogrammi
    // CosTheta [-1, 1]
    auto h_costheta_ptr = df.Histo1D({"h_costheta", ";Thrust cos(#theta);Events", 100, -1.0, 1.0}, "RP_thrustcostheta");
    // Phi [0, pi]
    auto h_phi_ptr = df.Histo1D({"h_phi", ";Thrust #phi;Events", 100, 0.0, TMath::Pi()}, "RP_thrustphi");

    TH1D *h_costheta = (TH1D*)h_costheta_ptr->Clone();
    TH1D *h_phi = (TH1D*)h_phi_ptr->Clone();

    // 3. Stile
    h_costheta->SetLineColor(kBlue + 1);
    h_costheta->SetLineWidth(2);
    h_costheta->SetFillColorAlpha(kBlue + 1, 0.3);

    h_phi->SetLineColor(kRed + 1);
    h_phi->SetLineWidth(2);
    h_phi->SetFillColorAlpha(kRed + 1, 0.3);

    // --- ASSI Y ---
    // Minimo a 0 e spazio sopra per la legenda
    h_costheta->SetMinimum(0); 
    //h_costheta->SetMaximum(h_costheta->GetMaximum() * 1.35); 

    h_phi->SetMinimum(0); 
    h_phi->SetMaximum(h_phi->GetMaximum() * 1.20); 

    // Setup dimensioni label (standard leggibili)
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

    // 4. Canvas
    TCanvas *c = new TCanvas("c", "Thrust Variables", 900, 800);
    // Divisione standard 1x2 (con margini di default tra i pad)
    c->Divide(1, 2); 

    // --- PAD 1 (CosTheta) ---
    c->cd(1);
    // Non forziamo margini strani, lasciamo quelli di default (o leggermente aggiustati per pulizia)
    //gPad->SetRightMargin(0.05); 
    //gPad->SetLeftMargin(0.12);
    //gPad->SetTopMargin(0.08); // Un po' di aria sopra
    
    h_costheta->Draw("HIST");

    // Legenda centrata in alto
    TLegend *leg1 = new TLegend(0.35, 0.75, 0.65, 0.92);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.06);
    leg1->SetTextAlign(22); 
    leg1->AddEntry(h_costheta, "Thrust cos(#theta)", "f");
    leg1->Draw();

    // --- PAD 2 (Phi) ---
    c->cd(2);
    //gPad->SetRightMargin(0.05);
    //gPad->SetLeftMargin(0.12);
    //gPad->SetTopMargin(0.08);
    //gPad->SetBottomMargin(0.12); // Spazio per l'etichetta dell'asse X
    
    h_phi->Draw("HIST");

    // Legenda centrata in alto
    TLegend *leg2 = new TLegend(0.35, 0.80, 0.65, 0.92);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.06);
    leg2->SetTextAlign(22);
    leg2->AddEntry(h_phi, "Thrust #phi", "f");
    leg2->Draw();

    // 5. Salvataggio
    gSystem->Exec(("mkdir -p " + outputDir).c_str());
    c->SaveAs(outputFile.c_str());
    
    std::cout << "Plot salvato in: " << outputFile << std::endl;
}
