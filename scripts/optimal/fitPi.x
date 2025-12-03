#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TString.h"
#include "TMatrixDSym.h"
#include "Fit/FitResult.h"

//================================================
// global functions
//================================================

// par[0]..par[7] : Pol7 bkg shape
// par[8]         : normalization (fraction f of background)
// par[9]         : N(+)
// par[10]        : N(-) 

double GlobalFitFunction(double *x, double *par) {
    
    double xx = x[0];
    double bkg_shape = 0;
    // Pol7: p0 + p1*x + ... + p7*x^7
    for(int i=0; i<=7; ++i) {
        bkg_shape += par[i] * pow(xx, i);
    }
    double BKG = par[8] * bkg_shape;

    // signal
    // w(+) = 1 + 1*(2x - 1) = 2x
    double w_plus = 2.0 * xx;
    
    // w(-) = 1 - 1*(2x - 1) = 2 - 2x
    double w_minus = 2.0 * (1.0 - xx);

    double SIG = par[9] * w_plus + par[10] * w_minus;
    double bin_width = par[11];
    double N = par[12];

    return N*bin_width*(BKG + SIG);
}

//================================================
// fitting function
//================================================
void Fit_Pi(TH1 *hist) {
        
    if (!hist) {
        cerr << "[ERROR] Fit_Pion_Global: hist is null\n";
        return;
    }

    // normalization
    TH1 *h = (TH1*)hist->Clone(Form("%s_clone", hist->GetName()));
    h->SetDirectory(0); 
    h->Sumw2();

    double N = h->Integral();
    double integral = h->Integral("width");
    double bin_width = h->GetBinWidth(1);
    

    // loading bkg_shape parameters
    TFile *inBKG = TFile::Open("bkg_parameters.root", "READ");
    if (!inBKG || inBKG->IsZombie()) {
        cerr << "[ERROR] Cannot open bkg_parameters.root. Run background fit";
        return;
    }
    
    TF1 *fBKG = (TF1*)inBKG->Get("bkg_shape");
    if (!fBKG) {
        cerr << "[ERROR] Function 'bkg_shape' not found in file.\n";
        inBKG->Close(); return;
    }

    // Global fit
    // 8 (Pol7) + 1 (Norm Bkg) + 1 (N+) + 1 (N-) = 11 params
    double xmin = 0.05;
    double xmax = 0.95;
    TF1 *f_tot = new TF1("gloabl_fit", GlobalFitFunction, xmin, xmax, 13);
    
    cout << "\n =========================\n" << endl;
	cout << "[INFO] Fitting " << hist->GetName() << " with " << N << " entries." << endl;

    // FIXING bkg parameters
    for(int i=0; i<=7; ++i) {
        double val = fBKG->GetParameter(i);
        //f_tot->SetParameter(i, val);
        f_tot->FixParameter(i, val); 
        f_tot->SetParName(i, Form("p%d_bkg", i));
    }
    inBKG->Close(); 

    //free parameters
    f_tot->SetParName(8,  "f_bkg"); // f , now means scale with respect to bkg fit (may have less events)
    f_tot->SetParName(9,  "N_{+}");   // N(+)
    f_tot->SetParName(10, "N_{-}");  // N(-)
    f_tot->SetParName(11,"bin_width");
    f_tot->SetParName(12,"entries");

    // initialization
    f_tot->SetParameter(8,  1); 
    f_tot->SetParameter(9,  N/2); 
    f_tot->SetParameter(10, N/2);
    f_tot->FixParameter(11,bin_width);
    f_tot->FixParameter(12,N);
    
    // positive and limited 
    //f_tot->SetParLimits(8,  0.0, 10.0);
    //f_tot->SetParLimits(9,  0.0, 10.0);
    //f_tot->SetParLimits(10, 0.0, 10.0);

    
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    
    TCanvas *c = new TCanvas(Form("c_%s", hist->GetName()), "Global Fit", 900, 600);
    c->cd();
    
    h->SetTitle("#tau #rightarrow #pi #nu_{#tau} channel;x_{#pi};Events");
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.7);
    h->SetLineColor(kBlack);
    h->SetLineWidth(1);
    h->SetMinimum(0.);
    h->Draw("HIST E1");
    
    f_tot->SetLineColor(kMagenta);
    f_tot->SetLineWidth(2);
    f_tot->Draw("SAME");

    // saving for cov matrix
    TFitResultPtr r = h->Fit(f_tot, "R S");

	

    // P = (N+ - N-) / (N+ + N-)
    double N_plus  = f_tot->GetParameter(9);
    double N_minus = f_tot->GetParameter(10);
    double P_val   = (N_plus - N_minus) / (N_plus + N_minus);

    // error propagated on P
    // dP/dN+ = 2*N- / (N+ + N-)^2
    // dP/dN- = -2*N+ / (N+ + N-)^2
    double den     = (N_plus + N_minus) * (N_plus + N_minus);
    double dP_dNplus  =  2.0 * N_minus / den;
    double dP_dNminus = -2.0 * N_plus  / den;
    
    // covariance matrix
    TMatrixDSym cov = r->GetCovarianceMatrix();
    double var_Nplus  = cov(9, 9);   // inidices match params
    double var_Nminus = cov(10, 10);
    double cov_NpNm   = cov(9, 10);
    
    
    double sigma_P2 = (dP_dNplus * dP_dNplus * var_Nplus) + 
                      (dP_dNminus * dP_dNminus * var_Nminus) + 
                      (2.0 * dP_dNplus * dP_dNminus * cov_NpNm);
    double P_err = sqrt(fabs(sigma_P2));

    // results
    cout << "\n ========================================= \n" << endl;
    cout << "[INFO] Fit Results for " << hist->GetName() << ":" << endl;
    cout << "   Norm Bkg (f) : " << f_tot->GetParameter(8) << endl;
    cout << "   N(+)         : " << N_plus << endl;
    cout << "   N(-)         : " << N_minus << endl;
    cout << "\n ========================================= \n" << endl;
    cout << "   Polarization : " << P_val << " +/- " << P_err << endl;
    double dev = fabs((fabs(P_val)-0.1468))/P_err;
    cout << "   P - P_th : " << setprecision(2) << dev << " std devs" << endl;
    cout << "\n ========================================= \n" << endl;
    
    // fit statistics
    double chi2 = f_tot->GetChisquare();  // Oppure r->Chi2()
    double ndf  = f_tot->GetNDF();        // Oppure r->Ndf()
    double prob = f_tot->GetProb();       // P-value

    cout << "\n ========================================= \n" << endl;
    cout << "   Chi2         : " << chi2 << endl;
    cout << "   NDF          : " << ndf << endl;
    cout << "   Chi2 / NDF   : " << chi2/ndf << endl;
    cout << "   Prob         : " << prob << endl;
    cout << "\n ========================================= \n" << endl;

    // Additional drawings
    TF1 *f_bkg = new TF1("bkg_draw", GlobalFitFunction, xmin, xmax, 13);
    f_bkg->SetParameters(f_tot->GetParameters());
    f_bkg->SetParameter(9, 0); // no N+
    f_bkg->SetParameter(10,0); // no N-
    f_bkg->SetLineColor(kGray+2);
    f_bkg->SetLineStyle(2);
    f_bkg->Draw("SAME");
    
    // Helicity +1 -1
    TF1 *H_plus = new TF1("Hplus_draw", GlobalFitFunction, xmin, xmax, 13);
    H_plus->SetParameters(f_tot->GetParameters());
    H_plus->SetParameter(8, 0); // no f
    H_plus->SetParameter(10,0); // no N-
    H_plus->SetLineColor(kBlue);
    H_plus->SetLineStyle(2);
    H_plus->Draw("SAME");
    TF1 *H_minus = new TF1("Hminus_draw", GlobalFitFunction, xmin, xmax, 13);
    H_minus->SetParameters(f_tot->GetParameters());
    H_minus->SetParameter(8, 0); // no f
    H_minus->SetParameter(9,0); // no N+
    H_minus->SetLineColor(kRed);
    H_minus->SetLineStyle(2);
    H_minus->Draw("SAME");
    
    // legend
    TLegend* leg = new TLegend(0.60,0.60,0.88,0.88);
    leg->AddEntry(h, "signal", "l");
    leg->AddEntry(f_tot, "global fit", "l");
    leg->AddEntry(f_bkg, "estimated background", "l");
    leg->AddEntry(H_plus, "H = +1", "l");
    leg->AddEntry(H_minus, "H = -1", "l");
    
    TString labelP = Form("P_{#tau} = %.4f #pm %.4f", P_val, P_err);
    leg->AddEntry((TObject*)0, labelP, "");
    leg->Draw();

    gPad->Update();

    // saving canvas
    const char* outdir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/plots/optimal/";
    c->SaveAs(TString(outdir) + "pion_fit.png");
    c->SaveAs(TString(outdir) + "pion_fit.pdf");

    // Pulizia
    delete f_bkg;
    delete H_plus;
    delete H_minus;
    delete h;
    delete c;
    delete f_tot;
}

//================================================
// MAIN MACRO
//================================================
void fitPi() {

    const char* infile = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/histmaker/optimal/p8_ee_Ztautau_ecm91.root";

    TFile* f = TFile::Open(infile, "READ");
    if (!f || f->IsZombie()) {
        cerr << "ERROR: cannot open file: " << infile << endl;
        return;
    }

    // Recupero istogramma pione (pi_sgn)
    TH1F* h_pi = (TH1F*)f->Get("pi_sgn");
	Fit_Pi(h_pi);
	
    //f->Close();
}
