#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>

void hadron_bkg() {
    // 1. Setup Paths
    std::string inputPath = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/treemaker/p8_ee_Ztautau_ecm91.root";
    const char* outDir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/RECO/additional/";

    // Ensure output directory exists
    gSystem->mkdir(outDir, true);

    // 2. Load Data 
    ROOT::EnableImplicitMT(); // Speed up reading
    ROOT::RDataFrame df("events", inputPath);

    // 3. Define Histogram Model
    // Range 0 to 47 GeV with bin width 0.5 -> (47-0)/0.5 = 94 bins
    auto histoModel = [](const char* name, const char* title) {
        return ROOT::RDF::TH1DModel(name, title, 94, 0, 47);
    };

    // =======================================================
    // 4. Create Histograms 
    // =======================================================
    
    // --- Pi Channel Variables ---
    auto h_ElAsPi  = df.Histo1D(histoModel("h_ElAsPi",  "Electron misid as #pi;Energy [GeV];Counts"), "ElAsPi_e");
    auto h_MuAsPi  = df.Histo1D(histoModel("h_MuAsPi",  "Muon misid as #pi;Energy [GeV];Counts"), "MuAsPi_e");
    auto h_RhoAsPi = df.Histo1D(histoModel("h_RhoAsPi", "#rho misid as #pi;Energy [GeV];Counts"), "RhoAsPi_e");
    auto h_A1AsPi  = df.Histo1D(histoModel("h_A1AsPi",  "a_{1} misid as #pi;Energy [GeV];Counts"), "A1AsPi_e");
    auto h_PiAsPi  = df.Histo1D(histoModel("h_PiAsPi",  "True #pi (Signal);Energy [GeV];Counts"), "PiAsPi_e");

    // --- Rho Channel Variables ---
    auto h_PiAsRho  = df.Histo1D(histoModel("h_PiAsRho",  "#pi misid as #rho;Energy [GeV];Counts"), "PiAsRho_e");
    auto h_RhoAsRho = df.Histo1D(histoModel("h_RhoAsRho", "True #rho (Signal);Energy [GeV];Counts"), "RhoAsRho_e");
    auto h_A1AsRho  = df.Histo1D(histoModel("h_A1AsRho",  "a_{1} misid as #rho;Energy [GeV];Counts"), "A1AsRho_e");
    auto h_AllRho   = df.Histo1D(histoModel("h_AllRho",   "Total Rho Candidates;Energy [GeV];Counts"), "AllRho_e");

    // =======================================================
    // 5. Styling and Plotting
    // =======================================================
    gStyle->SetOptStat(0); // Remove stat box
    
    // --- CANVAS 1: PION CHANNEL ---
    TCanvas *cPi = new TCanvas("cPi", "Pion Channel Backgrounds", 800, 600);
    TLegend *legPi = new TLegend(0.6, 0.6, 0.9, 0.9);
    
    // Helper lambda for styling
    auto styleHisto = [](TH1D* h, int color) {
        h->SetLineColor(color);
        h->SetLineWidth(2);
        h->SetFillColorAlpha(color, 0.35); // 35% Transparency
    };

    // Apply styles 
    styleHisto(h_PiAsPi.GetPtr(),  kBlack); // Signal
    styleHisto(h_ElAsPi.GetPtr(),  kRed);
    styleHisto(h_MuAsPi.GetPtr(),  kBlue);
    styleHisto(h_RhoAsPi.GetPtr(), kGreen+2);
    styleHisto(h_A1AsPi.GetPtr(),  kMagenta);

    // Draw
    h_PiAsPi->Draw("HIST"); // Draw signal first (or max) 
    // Ideally find max Y to set range, simplified here:
    h_PiAsPi->GetYaxis()->SetRangeUser(0, h_PiAsPi->GetMaximum()*1.2);
    
    h_ElAsPi->Draw("HIST SAME");
    h_MuAsPi->Draw("HIST SAME");
    h_RhoAsPi->Draw("HIST SAME");
    h_A1AsPi->Draw("HIST SAME");

    legPi->AddEntry(h_PiAsPi.GetPtr(), "True #pi", "f");
    legPi->AddEntry(h_ElAsPi.GetPtr(), "e #rightarrow #pi", "f");
    legPi->AddEntry(h_MuAsPi.GetPtr(), "#mu #rightarrow #pi", "f");
    legPi->AddEntry(h_RhoAsPi.GetPtr(), "#rho #rightarrow #pi", "f");
    legPi->AddEntry(h_A1AsPi.GetPtr(), "a_{1} #rightarrow #pi", "f");
    legPi->Draw();

    cPi->SaveAs(Form("%sPion_Backgrounds.pdf", outDir));

    // --- CANVAS 2: RHO CHANNEL ---
    TCanvas *cRho = new TCanvas("cRho", "Rho Channel Backgrounds", 800, 600);
    TLegend *legRho = new TLegend(0.6, 0.6, 0.9, 0.9);

    // Calculate "Rest of Background"
    // Rest = AllRho - (Pi + Rho + A1)
    TH1D *h_Rest = (TH1D*)h_AllRho->Clone("h_Rest");
    h_Rest->SetTitle("Other Backgrounds");
    h_Rest->Add(h_PiAsRho.GetPtr(), -1);
    h_Rest->Add(h_RhoAsRho.GetPtr(), -1);
    h_Rest->Add(h_A1AsRho.GetPtr(), -1);

    // Styling
    styleHisto(h_RhoAsRho.GetPtr(), kBlack);   // Signal
    styleHisto(h_PiAsRho.GetPtr(),  kRed);     // Backgrounds
    styleHisto(h_A1AsRho.GetPtr(),  kBlue);
    styleHisto(h_Rest,              kOrange+1); // Rest

    // Draw (Assumes Signal is largest, otherwise check max)
    h_RhoAsRho->Draw("HIST");
    h_RhoAsRho->GetYaxis()->SetRangeUser(0, h_RhoAsRho->GetMaximum()*1.2);

    h_PiAsRho->Draw("HIST SAME");
    h_A1AsRho->Draw("HIST SAME");
    h_Rest->Draw("HIST SAME");

    legRho->AddEntry(h_RhoAsRho.GetPtr(), "True #rho", "f");
    legRho->AddEntry(h_PiAsRho.GetPtr(), "#pi #rightarrow #rho", "f");
    legRho->AddEntry(h_A1AsRho.GetPtr(), "a_{1} #rightarrow #rho", "f");
    legRho->AddEntry(h_Rest, "Rest (Other Bkg)", "f");
    legRho->Draw();

    cRho->SaveAs(Form("%sRho_Backgrounds.pdf", outDir));

    // Cleanup memory
    delete legPi; delete cPi;
    delete legRho; delete cRho;
    delete h_Rest;
}