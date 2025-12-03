#include "TFractionFitter.h"
#include "TObjArray.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include <iostream>

void fitFraction() {

    const char* infile_data = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/treemaker/optimal/p8_ee_Ztautau_ecm91.root";
    const char* infile_templates = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/treemaker/optimal/bkg/templates.root";
    std::string treeName = "events";

    // recovering templates
    TFile *fTemp = TFile::Open(infile_templates);
    TH1D *h_plus  = (TH1D*)fTemp->Get("h_template_plus");
    TH1D *h_minus = (TH1D*)fTemp->Get("h_template_minus");
    h_plus->SetDirectory(0);
    h_minus->SetDirectory(0);
    fTemp->Close();

    // recover data and fill histogram
    ROOT::RDataFrame df(treeName, infile_data);
    auto h_data_ptr = df.Histo1D({"h_data", "Data", h_plus->GetNbinsX(), h_plus->GetXaxis()->GetXmin(), h_plus->GetXaxis()->GetXmax()}, "pi_sgn");
    TH1D *h_data = (TH1D*)h_data_ptr->Clone("h_data");
    
    
    h_data->Sumw2(); 

    // store TFraction histos
    TObjArray *mc = new TObjArray(2);
    mc->Add(h_plus);
    mc->Add(h_minus);
    
    TFractionFitter* fit = new TFractionFitter(h_data, mc);
    
    // FORSE SALTA QUI? BIN VUOTI GESTITI MALE?
	// nei primi bin il template -1 ha pesi negativi?
    fit->SetRangeX(h_data->FindBin(0.05), h_data->FindBin(0.95));

    // fitting
    // https://root.cern.ch/doc/master/classTFractionFitter.html
    Int_t status = fit->Fit();
    
    std::cout << "Fit Status: " << status << std::endl;

    if (status == 0) {
        
        double f_plus, f_minus, err_plus, err_minus;
        
        fit->GetResult(0, f_plus, err_plus);
        fit->GetResult(1, f_minus, err_minus);
        
        
        // P = (N+ - N-) / (N+ + N-) = (f+ - f-) / (f+ + f-)
        // TFractionFitter normalizes to f+ + f- = 1 ? TODO
        // P = f+ - f-
        
        // are fractions relative to the original templates?
        // we have to recover ->Integral()
        
        TH1F* resultPlot = (TH1F*)fit->GetPlot();
        
        // Calcolo manuale per sicurezza sulla normalizzazione
        double Int_plus_orig = h_plus->Integral();
        double Int_minus_orig = h_minus->Integral();
        
        // N_real = f_fit * (Integral_Data / Integral_Template)
        
        double P = (f_plus - f_minus) / (f_plus + f_minus);
        
        // error propagation (TODO add correlation from cov matrix)
        double err_P = sqrt(pow(2*f_minus/pow(f_plus+f_minus,2), 2)*pow(err_plus,2) + 
                            pow(-2*f_plus/pow(f_plus+f_minus,2), 2)*pow(err_minus,2));

        std::cout << "Fraction + : " << f_plus << " +/- " << err_plus << std::endl;
        std::cout << "Fraction - : " << f_minus << " +/- " << err_minus << std::endl;
        std::cout << "Polarization: " << P << " +/- " << err_P << std::endl;
        
        // Plot
        TCanvas *c = new TCanvas("c", "Fraction Fit", 800, 600);
        h_data->Draw("Ep");
        resultPlot->SetLineColor(kRed);
        resultPlot->Draw("same");
    }
}
