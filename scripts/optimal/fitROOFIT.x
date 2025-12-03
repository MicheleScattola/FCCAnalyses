#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooFit.h"
#include "ROOT/RDataFrame.hxx"
#include <iostream>

void fitROOFIT() {

    const char* infile_data = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/treemaker/optimal/p8_ee_Ztautau_ecm91.root";
    const char* infile_templates = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/treemaker/optimal/bkg/templates.root";
    std::string treeName = "events";

    // --- RECUPERO TEMPLATE ---
    TFile *fTemp = TFile::Open(infile_templates);
    TH1D *h_plus  = (TH1D*)fTemp->Get("h_template_plus");
    TH1D *h_minus = (TH1D*)fTemp->Get("h_template_minus");
    h_plus->SetDirectory(0);
    h_minus->SetDirectory(0);
    fTemp->Close();

    // --- HISTO DATA ---
    ROOT::RDataFrame df(treeName, infile_data);
    auto h_data_ptr = df.Histo1D(
        {"h_data", "Data", h_plus->GetNbinsX(),
         h_plus->GetXaxis()->GetXmin(), h_plus->GetXaxis()->GetXmax()},
        "pi_sgn"
    );
    TH1D *h_data = (TH1D*)h_data_ptr->Clone("h_data");
    h_data->Sumw2();

    // --- ROOFIT VARIABLE ---
    RooRealVar x("x","x",
                 h_plus->GetXaxis()->GetXmin(),
                 h_plus->GetXaxis()->GetXmax());

    // Limita range come in TFractionFitter
    x.setRange("fitRange", 0.05, 0.95);

    // --- ROO DATA HIST ---
    RooDataHist dh_data ("dh_data",  "data",  x, h_data);
    RooDataHist dh_plus ("dh_plus",  "plus",  x, h_plus);
    RooDataHist dh_minus("dh_minus", "minus", x, h_minus);

    // --- TEMPLATE PDFs ---
    RooHistPdf pdf_plus ("pdf_plus",  "pdf_plus",  x, dh_plus);
    RooHistPdf pdf_minus("pdf_minus", "pdf_minus", x, dh_minus);

    // --- FRAZIONE ---
    RooRealVar fplus("fplus","fraction plus", 0.5, 0.0, 1.0);

    // modello: f·H+ + (1−f)·H−
    RooAddPdf model("model","model",
                    RooArgList(pdf_plus, pdf_minus),
                    RooArgList(fplus));

    // --- FIT ---
    auto result = model.fitTo(dh_data,
                              RooFit::Range("fitRange"),
                              RooFit::Save(),
                              RooFit::PrintLevel(-1));

    double f_plus  = fplus.getVal();
    double err_plus = fplus.getError();
    double f_minus = 1.0 - f_plus;
    double err_minus = err_plus;   // completamente anticorrelati in un 2-template mixture

    // Polarizzazione
    double P = f_plus - f_minus;
    double err_P = sqrt(4 * err_plus * err_plus);

    std::cout << "Fraction + : " << f_plus  << " +/- " << err_plus  << std::endl;
    std::cout << "Fraction - : " << f_minus << " +/- " << err_minus << std::endl;
    std::cout << "Polarization: " << P << " +/- " << err_P << std::endl;

    // --- PLOT ---
    TCanvas *c = new TCanvas("c","c",800,600);
    //colors
    Int_t color_plus  = TColor::GetColorTransparent(kBlue, 0.5);
    Int_t color_minus = TColor::GetColorTransparent(kRed, 0.5);
    
    RooPlot* frame = x.frame();
    dh_data.plotOn(frame, RooFit::Name("data"));
    model.plotOn(frame, RooFit::Name("fit"), RooFit::LineColor(kBlack));
    model.plotOn(frame, RooFit::Components(pdf_plus),  RooFit::LineColor(kRed));
    model.plotOn(frame, RooFit::Components(pdf_minus), RooFit::LineColor(kBlue));
    frame->Draw();
}

