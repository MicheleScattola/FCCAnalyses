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
    double A_tau = par[0];
    double A_e = par[1]; // check universality
    
    double numerator = A_tau * (1 + costheta * costheta) + 2 * A_e * costheta;
    double denominator = 1 + costheta * costheta + 2 * A_tau * A_e * costheta;
    
    return -numerator / denominator;
  }
};

// Function to fit and plot a single TGraph
void plotChannel(const std::string &graph_name, const std::string &channel_label, 
                 const std::string &input_file, const std::string &output_pdf, 
                 const std::string &output_png) {
  
  gStyle->SetOptStat(0);
  
  // Open the file with the TGraph
  TFile *infile = TFile::Open(input_file.c_str(), "READ");
  if (!infile || infile->IsZombie()) {
    std::cerr << "[ERROR] Cannot open " << input_file << std::endl;
    return;
  }
  
  // Read the TGraph
  TGraphErrors *gr_Ptau = (TGraphErrors*)infile->Get(graph_name.c_str());
  if (!gr_Ptau) {
    std::cerr << "[ERROR] TGraph '" << graph_name << "' not found in file" << std::endl;
    infile->Close();
    return;
  }
  
  std::cout << "\n========================================" << std::endl;
  std::cout << "[INFO] Processing channel: " << channel_label << std::endl;
  std::cout << "[INFO] Loaded TGraph with " << gr_Ptau->GetN() << " points" << std::endl;
  
  // Restyle the graph
  gr_Ptau->SetTitle((channel_label + ";cos#theta;P_{#tau}").c_str());
  gr_Ptau->SetMarkerStyle(20);
  gr_Ptau->SetMarkerSize(1.0);
  gr_Ptau->SetMarkerColor(kBlack);
  gr_Ptau->SetLineColor(kBlack);
  gr_Ptau->SetLineWidth(2);
  gr_Ptau->GetYaxis()->SetRangeUser(-0.35, 0.05);
  gr_Ptau->GetXaxis()->SetRangeUser(-1, 1);
  
  // Create canvas
  TCanvas *c = new TCanvas("c_channel", ("Tau Polarization - " + channel_label).c_str(), 900, 700);
  c->SetGrid(1, 1);
  c->SetLeftMargin(0.14);
  
  // Create and fit with analytic function
  PtauFunctor functor;
  TF1 *f_ptau = new TF1("f_ptau", functor, -1.0, 1.0, 2);
  
  // Set parameter name and initial guess (start with SM value)
  f_ptau->SetParName(0, "A_{#tau}");
  f_ptau->SetParName(1, "A_{e}");
  f_ptau->SetParameter(0, SM_Atau_theory);
  f_ptau->SetParameter(1, SM_Atau_theory);
  
  // Fit to the data
  std::cout << ">>> Performing analytic fit for P_tau(cos theta)..." << std::endl;
  TFitResultPtr fitResult = gr_Ptau->Fit(f_ptau, "S");  // S = return TFitResultPtr
  
  // Extract results
  double A_tau_fit = f_ptau->GetParameter(0);
  double A_tau_err = f_ptau->GetParError(0);
  double A_e_fit = f_ptau->GetParameter(1);
  double A_e_err = f_ptau->GetParError(1);
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
  TF1 *f_theory = new TF1("f_theory", functor, -1.0, 1.0, 2);
  f_theory->SetParameter(0, SM_Atau_theory);
  f_theory->SetParameter(1, SM_Atau_theory);
  f_theory->SetLineColor(kBlue);
  f_theory->SetLineWidth(2);
  f_theory->SetLineStyle(2);
  f_theory->Draw("SAME");

  // Redraw data with error bars on top
  gr_Ptau->Draw("E1 P SAME");
  
  
  // Legend
  TLegend *leg = new TLegend(0.55, 0.55, 0.89, 0.89);
  leg->SetBorderSize(2);
  leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite);
  leg->SetTextSize(0.032);
  //leg->SetShadowColor(kGray + 2);
  leg->AddEntry(gr_Ptau, "Data", "ep");
  leg->AddEntry(f_ptau, Form("Fit: A_{#tau} = %.4f #pm %.4f", A_tau_fit, A_tau_err), "l");
  leg->AddEntry(f_ptau, Form("Fit: A_{e} = %.4f #pm %.4f", A_e_fit, A_e_err), "l");
  leg->AddEntry(f_theory, Form("SM: A_{e} #equiv A_{#tau} = %.4f", SM_Atau_theory), "l");
  leg->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %.2f", chi2_ndf), "");
  leg->Draw();
  
  /* Print results
  std::cout << std::string(60, '=') << std::endl;
  std::cout << ">>> " << channel_label << " FIT RESULTS <<<" << std::endl;
  std::cout << std::string(60, '=') << std::endl;
  std::cout << "\nFitted A_e (= A_tau): " << A_e_fit << " +/- " << A_e_err << std::endl;
  std::cout << "SM prediction A_e:    " << SM_Atau_theory << std::endl;
  std::cout << "Difference:           " << (A_e_fit - SM_Atau_theory) << std::endl;
  std::cout << "Relative deviation:   " << 100*(A_e_fit - SM_Atau_theory)/SM_Atau_theory << " %" << std::endl;
  std::cout << "\nFit Quality:" << std::endl;
  std::cout << "  χ²/ndf: " << chi2_ndf << std::endl;
  std::cout << "  χ²:     " << chi2 << std::endl;
  std::cout << "  ndf:    " << ndf << std::endl;
  std::cout << std::string(60, '=') << std::endl << std::endl;*/
  
  // Save the canvas
  c->SaveAs(output_pdf.c_str());
  c->SaveAs(output_png.c_str());
  std::cout << "[INFO] Canvas saved as " << output_pdf << " and " << output_png << std::endl;
  
  delete c;
  delete f_ptau;
  delete f_theory;
  delete leg;
  infile->Close();
}

void plot() {
  const std::string input_file = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/angular_analysis/universality.root";
  const std::string output_dir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/RECO/angular_analysis/";
  
  // Plot each channel
  plotChannel("el_theta", "Electron Channel", input_file, 
              output_dir + "universality_fit_el.pdf", output_dir + "universality_fit_el.png");
  
  plotChannel("mu_theta", "Muon Channel", input_file,
              output_dir + "universality_fit_mu.pdf", output_dir + "universality_fit_mu.png");
  
  plotChannel("pi_theta", "Pion Channel", input_file,
              output_dir + "universality_fit_pi.pdf", output_dir + "universality_fit_pi.png");
  
  plotChannel("rho_theta", "Rho Channel", input_file,
              output_dir + "universality_fit_rho.pdf", output_dir + "universality_fit_rho.png");
  
  plotChannel("combined_theta", "Combined Channel", input_file,
              output_dir + "universality_fit_combined.pdf", output_dir + "universality_fit_combined.png");
  
  std::cout << "\n>>> All plots completed successfully!" << std::endl;
}
