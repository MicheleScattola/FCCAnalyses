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
    // INTERNAL HELPERS (Hidden from header)
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
            double N = par[1];  // N of events, ideally fixed by histo integral
            double bw = par[2]; // bin_width
            
            double W = 1./3.*( (5 - 9*xx*xx + 4*xx*xx*xx) + P * (1 - 9*xx*xx + 8*xx*xx*xx) );
            double norm = (aux_norm(xmax, P) - aux_norm(xmin, P));
            
            if (norm == 0) return 0;
            return (N * bw * W) / norm;
        }
    };

    // =========================================================
    // IMPLEMENTATION 1: TEMPLATE FIT
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

        h_plus = (TH1D*)h_plus->Clone("mc_plus");
        h_minus = (TH1D*)h_minus->Clone("mc_minus");
        h_plus->SetDirectory(0); h_minus->SetDirectory(0);
        fTemp->Close();

        // 2. Get Data
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

        // 3. TFractionFitter
        TObjArray *mc = new TObjArray(2);
        mc->Add(h_plus); mc->Add(h_minus);
        TFractionFitter* fitter = new TFractionFitter(h_data, mc);
        fitter->Constrain(0, 0.0, 1.0);
        fitter->Constrain(1, 0.0, 1.0);

        Int_t status = fitter->Fit();
        if (status != 0) {
            std::cerr << "[Fitter] Template Fit failed." << std::endl;
            delete fitter; delete mc; delete h_data; delete h_plus; delete h_minus;
            return result;
        }

        // 4. Results
        double f_p, err_p, f_m, err_m;
        fitter->GetResult(0, f_p, err_p);
        fitter->GetResult(1, f_m, err_m);
        
        result.P_tau = f_p - f_m;
        result.P_err = 2.0 * err_p;
        result.f_plus = f_p;
        result.f_minus = f_m;
        result.success = true;

        std::cout << "[Fitter] " << plot_title << " (Template) | P_tau: " << result.P_tau << " +/- " << result.P_err << std::endl;

        // 5. Plotting (Simplified for brevity)
        gStyle->SetOptStat(0);
        TCanvas *c = new TCanvas("c", "Fit", 800, 600);
        TH1D* h_fit = (TH1D*)fitter->GetPlot();
        
        // Reconstruct components for drawing
        double scale_p = (f_p * h_data->Integral()) / h_plus->Integral();
        double scale_m = (f_m * h_data->Integral()) / h_minus->Integral();
        h_plus->Scale(scale_p); h_minus->Scale(scale_m);
        
        h_plus->SetLineColor(kBlue); h_plus->SetLineStyle(2); h_plus->SetFillColorALpha(kBlue, 0.1);
        h_minus->SetLineColor(kRed); h_minus->SetLineStyle(2); h_minus->SetFillColorAlpha(kRed, 0.1);
        h_fit->SetLineColor(kBlack); h_fit->SetLineWidth(2);
        h_data->SetMarkerStyle(20);

        // set y axis minimum to 0
        h_data->SetMinimum(0.);

        h_data->Draw("EP");
        h_fit->Draw("HIST SAME");
        h_plus->Draw("HIST SAME");
        h_minus->Draw("HIST SAME");

        TLegend *leg = new TLegend(0.6, 0.65, 0.88, 0.88);
        leg->AddEntry(h_data, "Data", "lp");
        leg->AddEntry(h_fit, "Global Fit", "l");
        leg->AddEntry(h_plus, "H=+1", "l");
        leg->AddEntry(h_minus, "H=-1", "l");
        leg->AddEntry((TObject*)0, Form("P = %.3f #pm %.3f", result.P_tau, result.P_err), "");
        leg->Draw();

        c->SaveAs((outdir + output_filename).c_str());
        
        delete c; delete fitter; delete mc; delete h_data; delete h_plus; delete h_minus;
        return result;
    }

    // =========================================================
    // IMPLEMENTATION 2: ANALYTIC FIT
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

        auto h_ptr = df.Histo1D({"h_data", (plot_title + ";" + x_axis_title + ";Events").c_str(), 20, 0.0, 1.0}, dataColName);
        TH1D *h_data = (TH1D*)h_ptr->Clone("h_analytic");
        h_data->SetDirectory(0);
        h_data->Sumw2();
        
        int N = h_data->Integral();
        if (N==0) return result;

        // 2. Analytic Function
        AnalyticFunctor functor(xmin, xmax);
        TF1 *f_fit = new TF1("f_ana", functor, xmin, xmax, 3);
        f_fit->SetParName(0, "P_tau");
        f_fit->SetParameter(0, -0.15);
        f_fit->FixParameter(1, N);
        f_fit->FixParameter(2, h_data->GetBinWidth(1));

        // 3. Fit
        gStyle->SetOptStat(0);
        h_data->Fit(f_fit, "R Q M E");

        result.P_tau = f_fit->GetParameter(0);
        result.P_err = f_fit->GetParError(0);
        result.success = true;

        std::cout << "[Fitter] " << plot_title << " (Analytic) | P_tau: " << result.P_tau << " +/- " << result.P_err << std::endl;

        // 4. Plotting
        TCanvas *c = new TCanvas("c_ana", "Analytic Fit", 800, 600);
        
        // Draw H+ / H- reference curves
        TF1 *h_p = new TF1("hp", functor, xmin, xmax, 3);
        h_p->SetParameters(1.0, N, h_data->GetBinWidth(1)); // P=+1
        h_p->SetLineColor(kBlue); h_p->SetLineStyle(2); h_p->SetFillColorAlpha(kBlue, 0.1);

        TF1 *h_m = new TF1("hm", functor, xmin, xmax, 3);
        h_m->SetParameters(-1.0, N, h_data->GetBinWidth(1)); // P=-1
        h_m->SetLineColor(kRed); h_m->SetLineStyle(2); h_m->SetFillColorAlpha(kRed, 0.1);

        h_data->SetMarkerStyle(20);
        h_data->Draw("EP");
        f_fit->Draw("SAME");
        h_p->Draw("SAME");
        h_m->Draw("SAME");

        TLegend *leg = new TLegend(0.2, 0.2, 0.45, 0.4);
        leg->AddEntry(h_data, "Data", "lp");
        leg->AddEntry(f_fit, "Fit", "l");
        leg->AddEntry(h_p, "H=+1", "l");
        leg->AddEntry(h_m, "H=-1", "l");
        leg->AddEntry((TObject*)0, Form("P = %.3f #pm %.3f", result.P_tau, result.P_err), "");
        leg->Draw();

        c->SaveAs((outdir + output_filename).c_str());

        delete c; delete f_fit; delete h_p; delete h_m; delete h_data;
        return result;
    }

} // namespace Fitter
