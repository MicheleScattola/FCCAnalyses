#include "FCCAnalyses/Fitter.h"

// ROOT Includes
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TFractionFitter.h"
#include "TObjArray.h"
#include "ROOT/RDataFrame.hxx"

namespace Fitter {

    // =========================================================
    // INTERNAL HELPERS 
    // =========================================================
    
    // Helper for Analytic Integral
    double aux_norm(double x, double P){
        return 1./3. * ( (5*x -3*x*x*x + x*x*x*x) + P * (x -3*x*x +2*x*x*x*x) );
    }

    // Functor for Analytic Fit
    struct AnalyticFunctor {
        const double xmin, xmax; 
        AnalyticFunctor(double min, double max) : xmin(min), xmax(max) {}

        double operator()(double *x, double *par) {
            double xx = x[0];
            double P = par[0];  // Polarization
            double N = par[1];  // N of events, fixed by histo integral
            double bw = par[2]; // bin_width
            
            double W = 1./3.*( (5 - 9*xx*xx + 4*xx*xx*xx) + P * (1 - 9*xx*xx + 8*xx*xx*xx) );
            double norm = (aux_norm(xmax, P) - aux_norm(xmin, P));
            
            if (norm == 0) return 0;
            return (N * bw * W) / norm;
        }
    };

    // Functor for Template Fit 
    // Data = Norm * [ ((1+P)/2)*H_plus + ((1-P)/2)*H_minus ]
    struct TemplateFunctor {
        TH1D *h_p, *h_m; 
        
        // Constructor takes pointers to the template histograms
        TemplateFunctor(TH1D* p, TH1D* m) : h_p(p), h_m(m) {}

        double operator()(double *x, double *par) {
            double xx = x[0];
            
            
            // We assume templates and data have identical binning
            int bin = h_p->FindBin(xx);
            
            double y_p = h_p->GetBinContent(bin);
            double y_m = h_m->GetBinContent(bin);

            // par[0] = Normalization
            // par[1] = Polarization P 
            double N = par[0];
            double P = par[1];

            // P = (N+ - N-) / (N+ + N-)
            // Coeffs: c+ = (1+P)/2, c- = (1-P)/2
            double c_plus = (1.0 + P) / 2.0;
            double c_minus = (1.0 - P) / 2.0;

            return N * (c_plus * y_p + c_minus * y_m);
        }
    };

    // NEW: Functor for Fraction Fit (Two Parameters N+, N-)
    // Data = N_plus * H_plus + N_minus * H_minus
    struct FractionFunctor {
        TH1D *h_p, *h_m; 
        FractionFunctor(TH1D* p, TH1D* m) : h_p(p), h_m(m) {}

        double operator()(double *x, double *par) {
            double xx = x[0];
            int bin = h_p->FindBin(xx);
            double y_p = h_p->GetBinContent(bin);
            double y_m = h_m->GetBinContent(bin);

            double N_plus = par[0];
            double N_minus = par[1];

            return N_plus * y_p + N_minus * y_m;
        }
    };

    // =========================================================
    // OLD METHOD: FRACTION FIT (Fit N+, N- then Error Prop)
    // =========================================================
    FitResult fit_fraction(const std::string& infile_data,
                           const std::string& infile_templates,
                           const std::string& outdir,
                           const std::string& output_filename, 
                           const std::string& treeName,
                           const std::string& dataColName,     
                           const std::string& name_plus,       
                           const std::string& name_minus,      
                           const std::string& plot_title,       
                           const std::string& x_axis_title      
                           ) 
    {
        FitResult result = {0,0,0,0,false, "Fraction"};

        std::cout << "[Fitter] Starting Fraction Fit (N+, N-) for " << plot_title << std::endl;

        // 1. Recover Templates
        TFile *fTemp = TFile::Open(infile_templates.c_str(), "READ");
        if (!fTemp || fTemp->IsZombie()) {
            std::cerr << "[Fitter] Error: Cannot open templates: " << infile_templates << std::endl;
            return result;
        }
        TH1D *h_plus  = (TH1D*)fTemp->Get(name_plus.c_str());
        TH1D *h_minus = (TH1D*)fTemp->Get(name_minus.c_str());
        
        if(!h_plus || !h_minus) { 
            std::cerr << "[Fitter] Template not found." << std::endl; return result; 
        }

        h_plus = (TH1D*)h_plus->Clone("mc_plus_frac");
        h_minus = (TH1D*)h_minus->Clone("mc_minus_frac");
        h_plus->SetDirectory(0); h_minus->SetDirectory(0);
        fTemp->Close();

        // Normalize templates to 1
        if (h_plus->Integral() > 0)  h_plus->Scale(1.0 / h_plus->Integral());
        if (h_minus->Integral() > 0) h_minus->Scale(1.0 / h_minus->Integral());

        // 2. Retrieve Data
        ROOT::EnableImplicitMT();
        ROOT::RDataFrame df(treeName, infile_data);
        if (!df.HasColumn(dataColName)) {
             std::cerr << "[Fitter] Data column missing: " << dataColName << std::endl; return result;
        }
        
        int nBins = h_plus->GetNbinsX();
        double xMin = h_plus->GetXaxis()->GetXmin();
        double xMax = h_plus->GetXaxis()->GetXmax();
        
        auto h_data_ptr = df.Histo1D({"h_data_frac", (plot_title + ";" + x_axis_title + ";Events").c_str(), nBins, xMin, xMax}, dataColName);
        TH1D *h_data = (TH1D*)h_data_ptr->Clone("data_frac");
        h_data->SetDirectory(0);
        h_data->Sumw2();

        // 3. Set up Fit Function
        FractionFunctor functor(h_plus, h_minus);
        TF1 *f_fit = new TF1("f_fraction", functor, xMin, xMax, 2);
        
        // Initial guesses: Split data 50/50
        double total_events = h_data->Integral();
        f_fit->SetParName(0, "N_plus");
        f_fit->SetParName(1, "N_minus");
        f_fit->SetParameter(0, total_events / 2.0);
        f_fit->SetParameter(1, total_events / 2.0);
        
        // 4. Perform Fit
        TFitResultPtr fitStatus = h_data->Fit(f_fit, "L S Q M E"); // S needed for Covariance Matrix access

        if ((Int_t)fitStatus != 0) {
            std::cerr << "[Fitter] Fraction Fit failed." << std::endl;
             return result;
        }

        // 5. Extract Results and Error Propagation
        double Np = f_fit->GetParameter(0);
        double Nm = f_fit->GetParameter(1);
        double eNp = f_fit->GetParError(0);
        double eNm = f_fit->GetParError(1);
        
        // Covariance
        TMatrixDSym cov = fitStatus->GetCovarianceMatrix();
        double covNpNm = cov(0,1);

        double N_tot = Np + Nm;
        
        // Polarization P = (Np - Nm) / (Np + Nm)
        double P_val = (Np - Nm) / N_tot;

        // Error Propagation
        // dP/dNp = 2*Nm / (Np+Nm)^2
        // dP/dNm = -2*Np / (Np+Nm)^2
        double dP_dNp = (2.0 * Nm) / (N_tot * N_tot);
        double dP_dNm = (-2.0 * Np) / (N_tot * N_tot);

        double var_P = (dP_dNp * dP_dNp * eNp * eNp) + 
                       (dP_dNm * dP_dNm * eNm * eNm) + 
                       (2.0 * dP_dNp * dP_dNm * covNpNm);

        result.P_tau = P_val;
        result.P_err = sqrt(var_P);
        result.f_plus  = Np / N_tot;
        result.f_minus = Nm / N_tot;
        result.success = true;

        std::cout << "[Fitter] " << plot_title << " (Fraction N+/N-) | P_tau: " << result.P_tau << " +/- " << result.P_err << std::endl;

        // 6. Plotting
        gStyle->SetOptStat(0);
        TCanvas *c = new TCanvas("c_frac", "Fraction Fit", 800, 600);
        
        // Scale templates for plotting using fit results
        TH1D* h_plus_plot = (TH1D*)h_plus->Clone("h_plus_plot_f");
        TH1D* h_minus_plot = (TH1D*)h_minus->Clone("h_minus_plot_f");
        
        h_plus_plot->Scale(Np);
        h_minus_plot->Scale(Nm);
        
        h_plus_plot->SetLineColor(kBlue); h_plus_plot->SetLineStyle(2); h_plus_plot->SetFillColorAlpha(kBlue, 0.1);
        h_minus_plot->SetLineColor(kRed); h_minus_plot->SetLineStyle(2); h_minus_plot->SetFillColorAlpha(kRed, 0.1);
        
        h_data->SetLineColor(kBlack);
        h_data->SetMarkerColor(kBlack);
        h_data->SetMarkerStyle(20);
        h_data->SetMarkerSize(0.8);
        h_data->SetMinimum(0.);

        TH1D *h_result_total = (TH1D*)h_plus_plot->Clone("h_res_total_f");
        h_result_total->Add(h_minus_plot); 
        
        h_result_total->SetLineColor(kGray+2);
        h_result_total->SetLineWidth(2);
        h_result_total->SetLineStyle(1);
        h_result_total->SetFillStyle(0);
        
        h_data->GetListOfFunctions()->Clear();
        
        h_data->Draw("E1 X0 P"); 
        h_result_total->Draw("HIST SAME");
        h_plus_plot->Draw("HIST SAME");
        h_minus_plot->Draw("HIST SAME");

        TLegend *leg = new TLegend(0.6, 0.6, 0.88, 0.9);
        leg->AddEntry(h_data, "Data", "lp");
        leg->AddEntry(h_plus_plot, "#it{H} = +1", "l");
        leg->AddEntry(h_minus_plot, "#it{H} = -1", "l");
        leg->AddEntry((TObject*)0, Form("#bf{#it{P} = %.4f #pm %.4f}", result.P_tau, result.P_err), "");
        leg->Draw();

        c->SaveAs((outdir + output_filename).c_str());
        
        delete c; delete f_fit; delete h_data; delete h_plus; delete h_minus; delete h_plus_plot; delete h_minus_plot;
        return result;
    }

    // =========================================================
    //TEMPLATE FIT (Using TF1 Linear Combo)
    // =========================================================
    FitResult fit(const std::string& infile_data,
                  const std::string& infile_templates,
                  const std::string& outdir,
                  const std::string& output_filename, 
                  const std::string& treeName,
                  const std::string& dataColName,     
                  const std::string& name_plus,       
                  const std::string& name_minus,      
                  const std::string& plot_title,      
                  const std::string& x_axis_title     
                  ) 
    {
        FitResult result = {0,0,0,0,false, "Template"};

        std::cout << "[Fitter] Starting Template Fit for " << plot_title << std::endl;

        //Recover Templates
        TFile *fTemp = TFile::Open(infile_templates.c_str(), "READ");
        if (!fTemp || fTemp->IsZombie()) {
            std::cerr << "[Fitter] Error: Cannot open templates: " << infile_templates << std::endl;
            return result;
        }
        TH1D *h_plus  = (TH1D*)fTemp->Get(name_plus.c_str());
        TH1D *h_minus = (TH1D*)fTemp->Get(name_minus.c_str());
        
        if(!h_plus || !h_minus) { 
            std::cerr << "[Fitter] Template not found." << std::endl; return result; 
        }

        h_plus = (TH1D*)h_plus->Clone("mc_plus");
        h_minus = (TH1D*)h_minus->Clone("mc_minus");
        h_plus->SetDirectory(0); h_minus->SetDirectory(0);
        fTemp->Close();

        // normalizing to 1 for safety (templates should already be PDFs)
        if (h_plus->Integral() > 0)  h_plus->Scale(1.0 / h_plus->Integral());
        if (h_minus->Integral() > 0) h_minus->Scale(1.0 / h_minus->Integral());

        // retrieve data
        ROOT::EnableImplicitMT();
        ROOT::RDataFrame df(treeName, infile_data);
        if (!df.HasColumn(dataColName)) {
             std::cerr << "[Fitter] Data column missing: " << dataColName << std::endl; return result;
        }
        
        int nBins = h_plus->GetNbinsX();
        double xMin = h_plus->GetXaxis()->GetXmin();
        double xMax = h_plus->GetXaxis()->GetXmax();
        
        auto h_data_ptr = df.Histo1D({"h_data", (plot_title + ";" + x_axis_title + ";Events").c_str(), nBins, xMin, xMax}, dataColName);
        TH1D *h_data = (TH1D*)h_data_ptr->Clone("data");
        h_data->SetDirectory(0);
        h_data->Sumw2();

        // linear fit with functor
        TemplateFunctor functor(h_plus, h_minus);
        TF1 *f_fit = new TF1("f_template", functor, xMin, xMax, 2);
        
        // Setup Parameters
        f_fit->SetParName(0, "Norm");
        f_fit->SetParName(1, "P_tau");
        f_fit->FixParameter(0, h_data->Integral());
        f_fit->SetParameter(1, -0.15); 
        
    
        // L = Log Likelihood
        // Q = Quiet
        // M = Improved Errors (alg from TMinuit)
        // E = Better errors estimation (Minos technique)
        TFitResultPtr fitStatus = h_data->Fit(f_fit, "L Q M E ");

        if ((Int_t)fitStatus != 0) {
            std::cerr << "[Fitter] Template Fit failed." << std::endl;
             // clean up...
             return result;
        }

        // results
        double P_val = f_fit->GetParameter(1);
        double P_err = f_fit->GetParError(1);
        double Norm  = f_fit->GetParameter(0);

        result.P_tau = P_val;
        result.P_err = P_err;
        // event fractions
        result.f_plus  = (1.0 + P_val) / 2.0; 
        result.f_minus = (1.0 - P_val) / 2.0;
        result.success = true;

        std::cout << "[Fitter] " << plot_title << " (Template TF1) | P_tau: " << result.P_tau << " +/- " << result.P_err << std::endl;

        // plotting
        gStyle->SetOptStat(0);
        TCanvas *c = new TCanvas("c", "Fit", 800, 600);
        
        // scale templates to data
        TH1D* h_plus_plot = (TH1D*)h_plus->Clone("h_plus_plot");
        TH1D* h_minus_plot = (TH1D*)h_minus->Clone("h_minus_plot");
        
        double scale_p = Norm * result.f_plus;
        double scale_m = Norm * result.f_minus;
        
        h_plus_plot->Scale(scale_p);
        h_minus_plot->Scale(scale_m);
        
        h_plus_plot->SetLineColor(kBlue); h_plus_plot->SetLineStyle(2); h_plus_plot->SetFillColorAlpha(kBlue, 0.1);
        h_minus_plot->SetLineColor(kRed); h_minus_plot->SetLineStyle(2); h_minus_plot->SetFillColorAlpha(kRed, 0.1);
        
        h_data->SetLineColor(kBlack);
        h_data->SetMarkerColor(kBlack);
        h_data->SetMarkerStyle(20);
        h_data->SetMarkerSize(0.8);
        h_data->SetMinimum(0.);

        TH1D *h_result_total = (TH1D*)h_plus_plot->Clone("h_res_total");
        h_result_total->Add(h_minus_plot); 
        
        h_result_total->SetLineColor(kGray+2);
        h_result_total->SetLineWidth(2);
        h_result_total->SetLineStyle(1);
        h_result_total->SetFillStyle(0);
        
        // remove fit function from histogram to avoid drawing it
        h_data->GetListOfFunctions()->Clear();
        
        h_data->Draw("E1 X0 P"); 
        h_result_total->Draw("HIST SAME");
        h_plus_plot->Draw("HIST SAME");
        h_minus_plot->Draw("HIST SAME");

        TLegend *leg = new TLegend(0.6, 0.6, 0.88, 0.9);
        leg->AddEntry(h_data, "Data", "lp");
        //leg->AddEntry(f_fit, "Global Fit", "l");
        leg->AddEntry(h_plus_plot, "#it{H} = +1", "l");
        leg->AddEntry(h_minus_plot, "#it{H} = -1", "l");
        leg->AddEntry((TObject*)0, Form("#bf{#it{P} = %.4f #pm %.4f}", result.P_tau, result.P_err), "");
        leg->Draw();

        c->SaveAs((outdir + output_filename).c_str());
        
        delete c; delete f_fit; delete h_data; delete h_plus; delete h_minus; delete h_plus_plot; delete h_minus_plot;
        return result;
    }

    // =========================================================
    // ANALYTIC FIT
    // =========================================================
    FitResult fit(const std::string& infile_data,
                  const std::string& outdir,
                  const std::string& output_filename, 
                  const std::string& treeName,
                  const std::string& dataColName,     
                  const std::string& plot_title,      
                  const std::string& x_axis_title,
                  double xmin,
                  double xmax
                  )
    {
        FitResult result = {0,0,0,0,false, "Analytic"};

        // 1. Get Data
        ROOT::EnableImplicitMT();
        ROOT::RDataFrame df(treeName, infile_data);
        if (!df.HasColumn(dataColName)) {
             std::cerr << "[Fitter] Data column missing: " << dataColName << std::endl; return result;
        }

        auto h_ptr = df.Histo1D({"h_data", (plot_title + ";" + x_axis_title + ";Events").c_str(), 40, 0.0, 1.0}, dataColName);
        TH1D *h_data = (TH1D*)h_ptr->Clone("h_analytic");
        h_data->SetDirectory(0);
        h_data->Sumw2();
        
        int N = h_data->Integral();
        if (N==0) return result;

        // analytic functor for denominator normalization (fitting range)
        AnalyticFunctor functor(xmin, xmax);
        TF1 *f_fit = new TF1("f_ana", functor, xmin, xmax, 3);
        f_fit->SetParName(0, "P_tau");
        f_fit->SetParameter(0, -0.15);
        f_fit->FixParameter(1, N);
        f_fit->FixParameter(2, h_data->GetBinWidth(1));

        // fit
        gStyle->SetOptStat(0);
        // R = use range
        h_data->Fit(f_fit, "L R Q M E");

        result.P_tau = f_fit->GetParameter(0);
        result.P_err = f_fit->GetParError(0);
        result.success = true;

        std::cout << "[Fitter] " << plot_title << " (Analytic) | P_tau: " << result.P_tau << " +/- " << result.P_err << std::endl;

        // plotting
        TCanvas *c = new TCanvas("c_ana", "Analytic Fit", 800, 600);
        
        // Draw H+ / H- reference curves
        TF1 *h_p = new TF1("hp", functor, xmin, xmax, 3);
        h_p->SetParameters(1.0, N, h_data->GetBinWidth(1)); // P=+1
        h_p->SetLineColor(kBlue); h_p->SetLineStyle(2); h_p->SetFillColorAlpha(kBlue, 0.1);

        TF1 *h_m = new TF1("hm", functor, xmin, xmax, 3);
        h_m->SetParameters(-1.0, N, h_data->GetBinWidth(1)); // P=-1
        h_m->SetLineColor(kRed); h_m->SetLineStyle(2); h_m->SetFillColorAlpha(kRed, 0.1);

        h_data->SetMarkerStyle(20);
        f_fit->SetLineColor(kMagenta);

        h_data->Draw("EP");
        f_fit->Draw("SAME");
        h_p->Draw("SAME");
        h_m->Draw("SAME");

        TLegend *leg = new TLegend(0.2, 0.2, 0.55, 0.4);
        leg->AddEntry(h_data, "Data", "lp");
        leg->AddEntry(f_fit, "Fit", "l");
        leg->AddEntry(h_p, "#it{H} = +1", "l");
        leg->AddEntry(h_m, "#it{H} = -1", "l");
        leg->AddEntry((TObject*)0, Form("#bf{#it{P} = %.4f #pm %.4f}", result.P_tau, result.P_err), "");
        leg->Draw();

        c->SaveAs((outdir + output_filename).c_str());

        delete c; delete f_fit; delete h_p; delete h_m; delete h_data;
        return result;
    }

} // namespace Fitter
