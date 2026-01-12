#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TFile.h"
#include <iostream>
#include <string>

// Helper function to handle plotting and saving for each particle
void create_and_save(ROOT::RDF::RNode df, 
                     TFile* fOut,
                     const std::string& var_sgn, 
                     const std::string& w_plus, 
                     const std::string& w_minus, 
                     const std::string& label, 
                     const std::string& suffix, 
                     const char* outdir,
                     int nBins, double xMin, double xMax) {

    // pointers to RResultPtr
    auto h_plus_ptr  = df.Histo1D({("h_plus_"+suffix).c_str(),  ("Template Helicity +1 (" + label + ");x;Events").c_str(), nBins, xMin, xMax}, var_sgn, w_plus);
    auto h_minus_ptr = df.Histo1D({("h_minus_"+suffix).c_str(), ("Template Helicity -1 (" + label + ");x;Events").c_str(), nBins, xMin, xMax}, var_sgn, w_minus);

    // Clone and detach from RDataFrame
    TH1D *h_plus  = (TH1D*)h_plus_ptr->Clone(("h_template_"+suffix+"_plus").c_str());
    TH1D *h_minus = (TH1D*)h_minus_ptr->Clone(("h_template_"+suffix+"_minus").c_str());

    // Normalization
    if (h_plus->Integral() > 0)  h_plus->Scale(1.0 / h_plus->Integral());
    if (h_minus->Integral() > 0) h_minus->Scale(1.0 / h_minus->Integral());

    cout << "[INFO] Created templates for " << label << ": " 
         << " Helicity +1 Integral = " << h_plus->Integral() 
         << ", Helicity -1 Integral = " << h_minus->Integral() << endl;

    // Styling
    h_plus->SetLineColor(kBlue);
    h_plus->SetLineWidth(2);
    h_plus->SetFillColorAlpha(kBlue, 0.3); 

    h_minus->SetLineColor(kRed);
    h_minus->SetLineWidth(2);
    h_minus->SetFillColorAlpha(kRed, 0.3);

    // Plotting
    TCanvas *c = new TCanvas(("c_" + suffix).c_str(), ("Polarization Templates " + label).c_str(), 900, 600);
    c->cd();

    // Determine Y range
    double max_y = std::max(h_plus->GetMaximum(), h_minus->GetMaximum());
    h_plus->SetMaximum(max_y * 1.25);
    h_plus->SetMinimum(0.);
    
    h_plus->SetTitle(("Polarization Templates " + label + ";x (re-weighted);Probability Density").c_str());

    h_plus->Draw("HIST");
    h_minus->Draw("HIST SAME");

    // Legend
    TLegend *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(h_plus,  "Helicity +1", "f");
    leg->AddEntry(h_minus, "Helicity -1", "f");
    leg->Draw();

    // save pdf
    c->SaveAs(TString(outdir) + "templates_" + suffix + ".pdf");

    // write to file
    fOut->cd();
    h_plus->Write();
    h_minus->Write();

    // Cleanup
    delete c;
    delete leg;
}

void templates() {
        
    // Settings
    const char* infile = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/treemaker/templates/p8_ee_Ztautau_ecm91.root";
    const char* outdir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/templates/";
    const char* outdir2 = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/treemaker/templates/";
    
    std::string treeName = "events"; 
    
    int nBins = 38;
    double xMin = 0.05;
    double xMax = 1.0; 

    // Style
    gStyle->SetOptStat(0); 

    // Open Dataframe
    ROOT::RDataFrame df(treeName, infile);

    // Open Output ROOT File (Single file for all histos)
    TString rootOutName = TString(outdir2) + "templates.root";
    TFile *fOut = new TFile(rootOutName, "RECREATE");

    std::cout << "Processing Pions..." << std::endl;
    create_and_save(df, fOut, "pi_sgn", "w_plus_pi", "w_minus_pi", "Pions", "pi", outdir, nBins, xMin, xMax);

    std::cout << "Processing Muons..." << std::endl;
    create_and_save(df, fOut, "mu_sgn", "w_plus_mu", "w_minus_mu", "Muons", "mu", outdir, nBins, xMin, xMax);

    std::cout << "Processing Electrons..." << std::endl;
    create_and_save(df, fOut, "el_sgn", "w_plus_el", "w_minus_el", "Electrons", "el", outdir, nBins, xMin, xMax);

    std::cout << "Processing Rho..." << std::endl;
    create_and_save(df, fOut, "rho_sgn", "w_plus_rho", "w_minus_rho", "Rho", "rho", outdir, 40, -1., 1.);


    // Close file
    fOut->Close();

    std::cout << "\n[INFO] All templates created and saved in: " << rootOutName << std::endl;
    std::cout << "[INFO] PDFs saved in: " << outdir << std::endl;
}
