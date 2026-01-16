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

    // Ensure output directory exists to avoid TPDF error
    gSystem->mkdir(outDir, true);

    // 2. Global Style Settings
    gStyle->SetOptTitle(kFALSE);      // Clean look: no title
    gStyle->SetOptStat(0);           // Clean look: no stats
    //gStyle->SetFillAlpha(0.35);      // Set transparency for Palette Fill Color (PFC)
    gStyle->SetPalette(kPastel);       

    // 3. Load Data
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df("events", inputPath);

    // Range 0 to 47 GeV with bin width 0.5 -> 94 bins
    auto hModel = [](const char* name, const char* title) {
        return ROOT::RDF::TH1DModel(name, title, 94, 0, 47);
    };

    // 4. Create Histograms (Lazy Evaluation)
    // The titles here will be used by gPad->BuildLegend()
    auto h_PiAsPi  = df.Histo1D(hModel("h_PiAsPi",  "True #pi"), "PiAsPi_e");
    auto h_ElAsPi  = df.Histo1D(hModel("h_ElAsPi",  "e #rightarrow #pi"), "ElAsPi_e");
    auto h_MuAsPi  = df.Histo1D(hModel("h_MuAsPi",  "#mu #rightarrow #pi"), "MuAsPi_e");
    auto h_RhoAsPi = df.Histo1D(hModel("h_RhoAsPi", "#rho #rightarrow #pi"), "RhoAsPi_e");
    auto h_A1AsPi  = df.Histo1D(hModel("h_A1AsPi",  "a_{1} #rightarrow #pi"), "A1AsPi_e");

    auto h_RhoAsRho = df.Histo1D(hModel("h_RhoAsRho", "True #rho"), "RhoAsRho_e");
    auto h_PiAsRho  = df.Histo1D(hModel("h_PiAsRho",  "#pi #rightarrow #rho"), "PiAsRho_e");
    auto h_A1AsRho  = df.Histo1D(hModel("h_A1AsRho",  "a_{1} #rightarrow #rho"), "A1AsRho_e");
    auto h_AllRho   = df.Histo1D(hModel("h_AllRho",   "Total Candidates"), "AllRho_e");

    // =======================================================
    // 5. Drawing - Pion Channel
    // =======================================================
    TCanvas *cPi = new TCanvas("cPi", "Pion Channel", 1000, 800);
    
    // Using HIST PLC PFC:
    // PLC picks Line Color from Palette
    // PFC picks Fill Color from Palette (transparency applied via gStyle)
    h_PiAsPi->GetXaxis()->SetTitle("Energy [GeV]");
    h_PiAsPi->GetYaxis()->SetTitle("Events")
    h_PiAsPi->GetYaxis()->SetRangeUser(0, h_ElAsPi->GetMaximum() * 1.1);

    h_PiAsPi->Draw("HIST PLC PFC"); 
    h_ElAsPi->Draw("HIST SAME PLC PFC");
    h_MuAsPi->Draw("HIST SAME PLC PFC");
    h_RhoAsPi->Draw("HIST SAME PLC PFC");
    h_A1AsPi->Draw("HIST SAME PLC PFC");

    cPi->Update();

    h_PiAsPi->SetFillColorAlpha(h_PiAsPi->GetFillColor(), 0.5);
    h_ElAsPi->SetFillColorAlpha(h_ElAsPi->GetFillColor(), 0.5);
    h_MuAsPi->SetFillColorAlpha(h_MuAsPi->GetFillColor(), 0.5);
    h_RhoAsPi->SetFillColorAlpha(h_RhoAsPi->GetFillColor(), 0.5);
    h_A1AsPi->SetFillColorAlpha(h_A1AsPi->GetFillColor(), 0.5);

    cPi->Modified();
    cPi->Update();

    gPad->BuildLegend(0.6, 0.65, 0.9, 0.9);
    cPi->SaveAs(Form("%sPion_Bkg.pdf", outDir));


    // =======================================================
    // 6. Drawing - Rho Channel
    // =======================================================
    gStyle->SetPalette(kLightTemperature); 
    TCanvas *cRho = new TCanvas("cRho", "Rho Channel", 1000, 800);

    // Calculate "Rest" background
    TH1D *h_Rest = (TH1D*)h_AllRho->Clone("h_Rest");
    h_Rest->SetTitle("Other bkg");
    h_Rest->Add(h_PiAsRho.GetPtr(), -1);
    h_Rest->Add(h_RhoAsRho.GetPtr(), -1);
    h_Rest->Add(h_A1AsRho.GetPtr(), -1);

    h_RhoAsRho->GetXaxis()->SetTitle("Energy [GeV]");
    h_RhoAsRho->GetYaxis()->SetTitle("Events / 0.5 GeV");
    h_RhoAsRho->GetYaxis()->SetRangeUser(0, h_RhoAsRho->GetMaximum() * 1.3);

    h_RhoAsRho->Draw("HIST PLC PFC");
    h_PiAsRho->Draw("HIST SAME PLC PFC");
    h_A1AsRho->Draw("HIST SAME PLC PFC");
    h_Rest->Draw("HIST SAME PLC PFC");

    cRho->Update();

    h_RhoAsRho->SetFillColorAlpha(h_RhoAsRho->GetFillColor(), 0.5);
    h_PiAsRho->SetFillColorAlpha(h_PiAsRho->GetFillColor(), 0.5);
    h_A1AsRho->SetFillColorAlpha(h_A1AsRho->GetFillColor(), 0.5);
    h_Rest->SetFillColorAlpha(h_Rest->GetFillColor(), 0.5);

    cRho->Modified();
    cRho->Update();

    gPad->BuildLegend(0.6, 0.65, 0.9, 0.9);
    cRho->SaveAs(Form("%sRho_Bkg.pdf", outDir));

    // 7. Cleanup
    delete cPi; delete cRho; delete h_Rest;
}
