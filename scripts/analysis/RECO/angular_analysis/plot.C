#include <iostream>
#include "TFile.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFitResultPtr.h"

// Constants from SM
const double SM_sin2thetaW = 0.2315;
const double gv_ga = 1 - 4 * SM_sin2thetaW;
const double SM_Atau_theory = 2 * gv_ga / (1 + gv_ga * gv_ga);

// Functor for analytic P_tau fit
// P_tau(cos theta) = -(A_tau * (1 + cos²θ) + 2*A_tau*cosθ) / (1 + cos²θ + 2*A_tau²*cosθ)
struct PtauFunctor {
  double operator()(double *x, double *par) const {
    double costheta = x[0];
    double A_tau = par[0];  // Parameter to fit
    
    double numerator = A_tau * (1 + costheta * costheta) + 2 * A_tau * costheta;
    double denominator = 1 + costheta * costheta + 2 * A_tau * A_tau * costheta;
    
    return -numerator / denominator;
  }
};

void plot() {
  gStyle->SetOptStat(0);
  //gStyle->SetOptFit(111);
  
  // Open the file with the TGraph
  TFile *infile = TFile::Open("/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/angular_analysis/universality.root", "READ");
  if (!infile || infile->IsZombie()) {
    std::cerr << "[ERROR] Cannot open universality.root" << std::endl;
    return;
  }
  
  // Read the TGraph
  TGraphErrors *gr_Ptau = (TGraphErrors*)infile->Get("graph_polarization");
  if (!gr_Ptau) {
    std::cerr << "[ERROR] TGraph 'graph_polarization' not found in file" << std::endl;
    infile->Close();
    return;
  }
  
  std::cout << "[INFO] Loaded TGraph with " << gr_Ptau->GetN() << " points" << std::endl;
  
  // Restyle the graph
  gr_Ptau->SetTitle("Tau Polarization vs cos(#theta);cos(#theta);P_{#tau}");
  gr_Ptau->SetMarkerStyle(20);
  gr_Ptau->SetMarkerSize(1.0);
  gr_Ptau->SetMarkerColor(kBlack);
  gr_Ptau->SetLineColor(kBlack);
  gr_Ptau->SetLineWidth(2);
  gr_Ptau->GetYaxis()->SetRangeUser(-0.35, 0.05);
  
  // Create canvas
  TCanvas *c = new TCanvas("c_universality", "Tau Polarization Universality", 900, 700);
  c->SetGrid(1, 1);
  c->SetLeftMargin(0.14);
  
  // Create and fit with analytic function
  PtauFunctor functor;
  TF1 *f_ptau = new TF1("f_ptau", functor, -1.0, 1.0, 1);
  
  // Set parameter name and initial guess (start with SM value)
  f_ptau->SetParName(0, "A_{e}");
  f_ptau->SetParameter(0, SM_Atau_theory);
  
  // Fit to the data
  std::cout << "\n>>> Performing analytic fit for P_tau(cos theta)..." << std::endl;
  TFitResultPtr fitResult = gr_Ptau->Fit(f_ptau, "S");  // S = return TFitResultPtr
  
  // Extract results
  double A_e_fit = f_ptau->GetParameter(0);
  double A_e_err = f_ptau->GetParError(0);
  double chi2 = f_ptau->GetChisquare();
  int ndf = f_ptau->GetNDF();
  double chi2_ndf = (ndf > 0) ? chi2 / ndf : -1;
  
  // Style the fit function
  f_ptau->SetLineColor(kRed);
  f_ptau->SetLineWidth(3);
  f_ptau->SetLineStyle(1);
  
  // Draw
  gr_Ptau->Draw("AP");
  f_ptau->Draw("SAME");
  
  // Add theory line
  TF1 *f_theory = new TF1("f_theory", functor, -1.0, 1.0, 1);
  f_theory->SetParameter(0, SM_Atau_theory);
  f_theory->SetLineColor(kBlue);
  f_theory->SetLineWidth(2);
  f_theory->SetLineStyle(2);
  f_theory->Draw("SAME");
  
  // Legend
  TLegend *leg = new TLegend(0.58, 0.60, 0.90, 0.88);
  leg->SetBorderSize(1);
  leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite);
  leg->AddEntry(gr_Ptau, "Combined data", "ep");
  leg->AddEntry(f_ptau, Form("Fit: A_{e} = %.4f #pm %.4f", A_e_fit, A_e_err), "l");
  leg->AddEntry(f_theory, Form("SM: A_{e} = %.4f", SM_Atau_theory), "l");
  leg->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %.2f", chi2_ndf), "");
  leg->Draw();
  
  // Print results
  std::cout << "\n" << std::string(60, '=') << std::endl;
  std::cout << ">>> UNIVERSALITY FIT RESULTS <<<" << std::endl;
  std::cout << std::string(60, '=') << std::endl;
  std::cout << "\nFitted A_e (= A_tau): " << A_e_fit << " +/- " << A_e_err << std::endl;
  std::cout << "SM prediction A_e:    " << SM_Atau_theory << std::endl;
  std::cout << "Difference:           " << (A_e_fit - SM_Atau_theory) << std::endl;
  std::cout << "Relative deviation:   " << 100*(A_e_fit - SM_Atau_theory)/SM_Atau_theory << " %" << std::endl;
  std::cout << "\nFit Quality:" << std::endl;
  std::cout << "  χ²/ndf: " << chi2_ndf << std::endl;
  std::cout << "  χ²:     " << chi2 << std::endl;
  std::cout << "  ndf:    " << ndf << std::endl;
  std::cout << std::string(60, '=') << std::endl << std::endl;
  
  // Save the canvas
  c->SaveAs("/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/RECO/universality_fit.pdf");
  c->SaveAs("/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/RECO/universality_fit.png");
  std::cout << "[INFO] Canvas saved as universality_fit.pdf and universality_fit.png" << std::endl;
  
  // Optionally save results to ROOT file
  TFile *outfile = TFile::Open("/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/angular_analysis/universality_results.root", "RECREATE");
  gr_Ptau->Write("gr_Ptau");
  f_ptau->Write("f_ptau_fit");
  f_theory->Write("f_ptau_theory");
  outfile->Close();
  std::cout << "[INFO] Results saved to universality_results.root" << std::endl;
  
  infile->Close();
}
