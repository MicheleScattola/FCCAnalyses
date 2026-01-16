#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TH1D.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <vector>

void hadron_bkg_markers() {
    // 1. Setup Paths
    std::string inputPath = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/treemaker/p8_ee_Ztautau_ecm91.root";
    const char* outDir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/RECO/additional/";

    // Ensure output directory exists
    gSystem->mkdir(outDir, true);

    // 2. Global Style Settings (From your snippet)
    gStyle->SetOptTitle(kFALSE); // No title on top
    gStyle->SetOptStat(0);       // No stat box
    gStyle->SetPalette(kBird);   // Default rainbow-like palette, works well with PLC

    // 3. Load Data
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df("events", inputPath);

    // 4. Define Histograms
    // Note: The Title string here becomes the Legend entry automatically!
    auto hModel = [](const char* name, const char* title) {
        return ROOT::RDF::TH1DModel(name, title, 94, 0, 47);
    };

    // --- Pi Channel ---
    auto h_PiAsPi  = df.Histo1D(hModel("h_PiAsPi",  "True #pi (Signal)"), "PiAsPi_e");
    auto h_ElAsPi  = df.Histo1D(hModel("h_ElAsPi",  "e #rightarrow #pi"), "ElAsPi_e");
    auto h_MuAsPi  = df.Histo1D(hModel("h_MuAsPi",  "#mu #rightarrow #pi"), "MuAsPi_e");
    auto h_RhoAsPi = df.Histo1D(hModel("h_RhoAsPi", "#rho #rightarrow #pi"), "RhoAsPi_e");
    auto h_A1AsPi  = df.Histo1D(hModel("h_A1AsPi",  "a_{1} #rightarrow #pi"), "A1AsPi_e");

    // --- Rho Channel ---
    auto h_RhoAsRho = df.Histo1D(hModel("h_RhoAsRho", "True #rho (Signal)"), "RhoAsRho_e");
    auto h_PiAsRho  = df.Histo1D(hModel("h_PiAsRho",  "#pi #rightarrow #rho"), "PiAsRho_e");
    auto h_A1AsRho  = df.Histo1D(hModel("h_A1AsRho",  "a_{1} #rightarrow #rho"), "A1AsRho_e");
    auto h_AllRho   = df.Histo1D(hModel("h_AllRho",   "Total Candidates"), "AllRho_e");

    // =======================================================
    // 5. Plotting Logic
    // =======================================================

    // Helper to apply marker styles
    auto setStyle = [](TH1D* h, int markerStyle) {
        h->SetMarkerStyle(markerStyle);
        h->SetMarkerSize(1.2); // Make them slightly visible
    };

    // --- CANVAS 1: PION CHANNEL ---
    TCanvas *cPi = new TCanvas("cPi", "Pion Channel", 1000, 800);
    
    // Apply Styles
    setStyle(h_PiAsPi.GetPtr(),  kFullCircle);
    setStyle(h_ElAsPi.GetPtr(),  kFullSquare);
    setStyle(h_MuAsPi.GetPtr(),  kFullTriangleUp);
    setStyle(h_RhoAsPi.GetPtr(), kFullTriangleDown);
    setStyle(h_A1AsPi.GetPtr(),  kOpenCircle);

    // Drawing with "PLC PMC" (Palette Line Color, Palette Marker Color)
    h_ElAsPi->Draw("P PLC PMC"); 
    h_PiAsPi->Draw("SAME P PLC PMC");
    h_MuAsPi->Draw("SAME P PLC PMC");
    h_RhoAsPi->Draw("SAME P PLC PMC");
    h_A1AsPi->Draw("SAME P PLC PMC");

    // Auto-build legend based on titles
    gPad->BuildLegend(0.6, 0.7, 0.9, 0.9); 

    cPi->SaveAs(Form("%sPion_Backgrounds_Markers.pdf", outDir));


    // --- CANVAS 2: RHO CHANNEL ---
    TCanvas *cRho = new TCanvas("cRho", "Rho Channel", 1000, 800);

    // Calculate "Rest"
    TH1D *h_Rest = (TH1D*)h_AllRho->Clone("h_Rest");
    h_Rest->SetTitle("Rest (Other Bkg)");
    h_Rest->Add(h_PiAsRho.GetPtr(), -1);
    h_Rest->Add(h_RhoAsRho.GetPtr(), -1);
    h_Rest->Add(h_A1AsRho.GetPtr(), -1);

    // Apply Styles
    setStyle(h_RhoAsRho.GetPtr(), kFullCircle);
    setStyle(h_PiAsRho.GetPtr(),  kFullSquare);
    setStyle(h_A1AsRho.GetPtr(),  kFullTriangleUp);
    setStyle(h_Rest,              kFullTriangleDown);

    // Draw
    h_RhoAsRho->GetYaxis()->SetRangeUser(0, h_RhoAsRho->GetMaximum() * 1.3);
    h_RhoAsRho->Draw("P PLC PMC");
    h_PiAsRho->Draw("SAME P PLC PMC");
    h_A1AsRho->Draw("SAME P PLC PMC");
    h_Rest->Draw("SAME P PLC PMC");

    // Auto-build legend
    gPad->BuildLegend(0.6, 0.7, 0.9, 0.9);

    cRho->SaveAs(Form("%sRho_Backgrounds_Markers.pdf", outDir));

    // Cleanup
    delete cPi;
    delete cRho;
    delete h_Rest;
}
