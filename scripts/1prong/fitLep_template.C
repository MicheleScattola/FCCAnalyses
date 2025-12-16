#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TString.h"
#include "TLegend.h"
#include "TMatrixDSym.h"
#include "Fit/FitResult.h"
#include "ROOT/RDataFrame.hxx"

// =============================================================================
// functor class to pass histograms
struct TemplateFitFunctor {
    const TH1D* hP; // Template Helicity +1
    const TH1D* hM; // Template Helicity -1

    // constructor
    TemplateFitFunctor(const TH1D* h_plus, const TH1D* h_minus) : hP(h_plus), hM(h_minus) {}


    double operator()(double *x, double *par) {
        double xx = x[0];

        double val_plus  = hP->Interpolate(xx); 
        double val_minus = hM->Interpolate(xx);
		// N(+) * h(+) + N(-) * h(-)
        return par[0] * val_plus + par[1] * val_minus;
    }
};

// =============================================================================

void fitLep_template() {

    const char* infile_data = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/treemaker/p8_ee_Ztautau_ecm91.root";
    
    const char* infile_templates = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/treemaker/bkg/templates.root";
    
    const char* outdir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/";
    std::string treeName = "events";

    // recover template
    TFile *fTemp = TFile::Open(infile_templates, "READ");
    
    TH1D *h_plus  = (TH1D*)fTemp->Get("h_template_el_plus");
    TH1D *h_minus = (TH1D*)fTemp->Get("h_template_el_minus");
    
    h_plus->SetDirectory(0);
    h_minus->SetDirectory(0);
    fTemp->Close();

    // create histogram of signal data
    int nBins = h_plus->GetNbinsX();
    double xMin = h_plus->GetXaxis()->GetXmin();
    double xMax = h_plus->GetXaxis()->GetXmax();

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df(treeName, infile_data);
    
    auto h_data_ptr = df.Histo1D({"h_data", "Fit Polarization;x_{e};Events", nBins, xMin, xMax}, "el_sgn");
    TH1D *h_data = (TH1D*)h_data_ptr->Clone("h_data_final");

    // fit using functor
    TemplateFitFunctor fitFunctor(h_plus, h_minus);
    
    TF1 *f_fit = new TF1("f_pol_fit", fitFunctor, 0.05, 1.0, 2);
    f_fit->SetParName(0, "N_plus");
    f_fit->SetParName(1, "N_minus");

    // initialization
    double integral_tot = h_data->Integral();
    f_fit->SetParameter(0, integral_tot * 0.4); 
    f_fit->SetParameter(1, integral_tot * 0.6); 
    f_fit->SetParLimits(0, 0.0, integral_tot * 2);
    f_fit->SetParLimits(1, 0.0, integral_tot * 2);

    // plot
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111); 

    TCanvas *c = new TCanvas("c_fit", "Polarization Fit", 800, 600);
    c->cd();
    TGaxis::SetMaxDigits(3);
	// S saves data for cov matrix 
    TFitResultPtr r = h_data->Fit(f_fit, "L S R"); 

    double N_plus  = f_fit->GetParameter(0);
    double N_minus = f_fit->GetParameter(1);
    
    // P = (N+ - N-) / (N+ + N-)
    double P_val = (N_plus - N_minus) / (N_plus + N_minus);
    
    // error propagation with cov matrix
    TMatrixDSym cov = r->GetCovarianceMatrix();
    double var_Np = cov(0, 0);
    double var_Nm = cov(1, 1);
    double cov_NpNm = cov(0, 1);
    
    double S = N_plus + N_minus; // sum S
    double D = N_plus - N_minus; // diff D
    
    // dP/dN+ = 2*N- / S^2
    // dP/dN- = -2*N+ / S^2
    double dP_dNp =  2.0 * N_minus / (S*S);
    double dP_dNm = -2.0 * N_plus  / (S*S);
    
    double sigma2_P = (dP_dNp*dP_dNp * var_Np) + 
                      (dP_dNm*dP_dNm * var_Nm) + 
                      (2.0 * dP_dNp * dP_dNm * cov_NpNm);
    
    double P_err = sqrt(sigma2_P);

    std::cout << "\n=============================================" << std::endl;
    std::cout << " FIT RESULTS " << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << " N(+) : " << N_plus << " +/- " << f_fit->GetParError(0) << std::endl;
    std::cout << " N(-) : " << N_minus << " +/- " << f_fit->GetParError(1) << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << " P_tau: " << P_val << " +/- " << P_err << std::endl;
    std::cout << "=============================================\n" << std::endl;
    
    // H(+)
    TH1D *h_result_plus = (TH1D*)h_plus->Clone("h_res_plus");
    //h_result_plus->Reset(); 
    //h_result_plus->Add(h_plus, N_plus);
    h_result_plus->Scale(N_plus);
    
    h_result_plus->SetLineColor(kBlue);
    h_result_plus->SetLineWidth(2);
    h_result_plus->SetLineStyle(2); 
    h_result_plus->SetFillColorAlpha(kBlue, 0.05); 

    // H(-)
    TH1D *h_result_minus = (TH1D*)h_minus->Clone("h_res_minus");
    //h_result_minus->Reset();
    //h_result_minus->Add(h_minus, N_minus); 
    h_result_minus->Scale(N_minus);
    
    h_result_minus->SetLineColor(kRed);
    h_result_minus->SetLineWidth(2);
    h_result_minus->SetLineStyle(2); 
    h_result_minus->SetFillColorAlpha(kRed, 0.05); 

    // Tot fit = sum
    TH1D *h_result_total = (TH1D*)h_result_plus->Clone("h_res_total");
    h_result_total->Add(h_result_minus); 
    
    h_result_total->SetLineColor(kBlack);
    h_result_total->SetLineWidth(2);
    h_result_total->SetLineStyle(1); 
    h_result_total->SetFillStyle(0); 

    h_data->SetLineColor(kBlack);
    h_data->SetMarkerColor(kBlack);
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.8);
    h_data->SetTitle("Tau Polarization Fit;x_{e};Events");
    h_data->SetMinimum(0.);
    

    double max_data = h_data->GetMaximum();
    double max_fit  = h_result_total->GetMaximum();
    h_data->SetMaximum(std::max(max_data, max_fit) * 1.2);


    h_data->Draw("EP HIST");             
    h_result_total->Draw("HIST SAME"); 
    h_result_plus->Draw("HIST SAME"); 
    h_result_minus->Draw("HIST SAME"); 
    
    // f_fit->Draw()

    // Legend
    //TLegend *leg = new TLegend(0.55, 0.6, 0.88, 0.88);
    TLegend *leg = new TLegend(0.65, 0.65, 0.88, 0.73);
    //leg->SetBorderSize(0);
    //leg->SetFillStyle(0); 
    //leg->AddEntry(h_data, "Data", "lep");
    //leg->AddEntry(h_result_total, "Global Fit (Total)", "l");
    //leg->AddEntry(h_result_plus,  "Helicity +1", "lf");
    //leg->AddEntry(h_result_minus, "Helicity -1", "lf");
    leg->AddEntry((TObject*)0, Form("P_{#tau} = %.4f #pm %.4f", P_val, P_err), "");
    
    leg->Draw();

    // save
    //c->SaveAs(TString(outdir) + "polarization_fit_result_histo.png");
    c->SaveAs(TString(outdir) + "el_template_fit.pdf");

    delete h_result_plus;
    delete h_result_minus;
    delete h_result_total;
}
