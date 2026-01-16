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
    //gStyle->SetPalette(kRedBlue);       
    TGaxis::SetMaxDigits(3);
    gStyle->SetPadRightMargin(0.05);

    // 3. Load Data
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df("events", inputPath);

    // Range 0 to 47 GeV with bin width 0.5 -> 94 bins
    auto hModel = [](const char* name, const char* title) {
        return ROOT::RDF::TH1DModel(name, title, 94, 0, 47);
    };
    // Helper lambda for styling
    auto styleHisto = [](TH1D* h, int color) {
        h->SetLineColor(color);
        h->SetLineWidth(2);
        h->SetFillColorAlpha(color, 0.4); 
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
     
    // set styles
    h_PiAsPi->GetXaxis()->SetTitle("Energy [GeV]");
    h_PiAsPi->GetYaxis()->SetTitle("Events");
    h_PiAsPi->GetYaxis()->SetRangeUser(0, h_ElAsPi->GetMaximum() * 1.1);

    styleHisto(h_PiAsPi.GetPtr(),  kGreen); 
    styleHisto(h_ElAsPi.GetPtr(),  kRed);
    styleHisto(h_MuAsPi.GetPtr(),  kMagenta);
    styleHisto(h_RhoAsPi.GetPtr(), kBlue);
    styleHisto(h_A1AsPi.GetPtr(),  KAzure);

    h_PiAsPi->Draw("HIST"); 
    h_ElAsPi->Draw("HIST SAME");
    h_MuAsPi->Draw("HIST SAME");
    h_RhoAsPi->Draw("HIST SAME");
    h_A1AsPi->Draw("HIST SAME");
    
    gPad->BuildLegend();
    gPad->RedrawAxis();

    cPi->SaveAs(Form("%sPion_Bkg.pdf", outDir));


    // =======================================================
    // 6. Drawing - Rho Channel
    // =======================================================
    //gStyle->SetPalette(kCool); 
    TCanvas *cRho = new TCanvas("cRho", "Rho Channel", 1000, 800);
   
    // Calculate "Rest" background
    TH1D *h_Rest = (TH1D*)h_AllRho->Clone("h_Rest");
    h_Rest->SetTitle("Other bkg");
    h_Rest->Add(h_PiAsRho.GetPtr(), -1);
    h_Rest->Add(h_RhoAsRho.GetPtr(), -1);
    h_Rest->Add(h_A1AsRho.GetPtr(), -1);

    // set styles
    h_RhoAsRho->GetXaxis()->SetTitle("Energy [GeV]");
    h_RhoAsRho->GetYaxis()->SetTitle("Events");
    h_RhoAsRho->GetYaxis()->SetRangeUser(0, h_RhoAsRho->GetMaximum() * 1.3);

    styleHisto(h_RhoAsRho.GetPtr(), kBlue);   
    styleHisto(h_PiAsRho.GetPtr(),  kGreen);     
    styleHisto(h_A1AsRho.GetPtr(),  kAzure);
    styleHisto(h_Rest.GetPtr(),     kGray+2);

    h_RhoAsRho->Draw("HIST");
    
    h_A1AsRho->Draw("HIST SAME");
    h_Rest->Draw("HIST SAME");
    h_PiAsRho->Draw("HIST SAME");
    
    gPad->BuildLegend();
    gPad->RedrawAxis();
    
    cRho->SaveAs(Form("%sRho_Bkg.pdf", outDir));

    // 7. Cleanup
    delete cPi; delete cRho; delete h_Rest;
}
