#include <iostream>
#include "ROOT/RDataFrame.hxx"
//================================================
double aux_norm(double x, double P){
	return 1./3. * ( (5*x -3*x*x*x + x*x*x*x) + P * (x -3*x*x +2*x*x*x*x) );
}
// =============================================================================
// functor class to pass fitting range normalization
struct TemplateFitFunctor {
    const double xmin; 
    const double xmax; 

    // constructor
    TemplateFitFunctor(const double &x_min, const double &x_max) : xmin(x_min), xmax(x_max) {}


    double operator()(double *x, double *par) {

        double W =  par[1]*par[2]*1./3.*( (5-9*x[0]*x[0]+4*x[0]*x[0]*x[0]) + par[0] * (1-9*x[0]*x[0]+8*x[0]*x[0]*x[0]) );
		double norm = ( aux_norm(xmax,par[0]) - aux_norm(xmin,par[0]) );
        return W/norm;
    }
};	

void Fit_1prong(TH1D *hist, const char *outname, const char *title) {
		
	if (!hist) {
		std::cerr << "[ERROR] Fit_1prong: hist is null\n";
		return;
	}

	// cloning and normalizing histogram
	TH1D *h = (TH1D*)hist->Clone(Form("%s_clone", hist->GetName()));
	h->Sumw2();
	
	int N = h->Integral();
	cout << "\n =========================\n" << endl;
	cout << "[INFO] Fitting " << hist->GetName() << " with " << N << " entries." << endl;
	
	double bin_width = h->GetBinWidth(1);
	
	double xmin = h->GetXaxis()->GetXmin();
	double xmax = h->GetXaxis()->GetXmax();
	
	// fit using functor
    TemplateFitFunctor fitFunctor(xmin, xmax);
    
    TF1 *f = new TF1("fit", fitFunctor, xmin, xmax, 3);
	cout << "Fitting between " << xmin << " and " << xmax << endl;     
	// initialize values
	f->SetParName(0, "P_{#tau}");
	f->SetParameter(0, -0.1); 
	f->FixParameter(1, N);
	f->FixParameter(2, bin_width);
	
	double Perr = 0.0;
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1111);
	
	
	// plot
	TCanvas *c = new TCanvas(Form("canvas_%s", hist->GetName()), "P fit", 900,600);
	c->cd();
	TGaxis::SetMaxDigits(3);
	// Fit:
	
	h->SetTitle(title);
	h->SetMarkerStyle(20);
	h->SetMarkerSize(0.7);
	h->SetLineColor(kBlack);
	h->SetLineWidth(2);
	      
	f->SetLineColor(kMagenta);
	f->SetLineWidth(2);
	
	h->SetStats(1);
	h->Draw("HIST E");
	h->Fit(f, "RQM SAME");
	f->Draw("SAME");
	
	// results
	double P = f->GetParameter(0);
	Perr = f->GetParError(0);

	std::cout << "[INFO] Fit obtained for " << hist->GetName() << ":" << endl;
	std::cout << "****** Polarization P = " << P << " +/- " << Perr << endl;
	std::cout << "PDF normalized at : " << aux_norm(xmax,P) - aux_norm(xmin,P) << endl;
	cout << "\n =========================\n" << endl;
	
	// Additional drawings
    
    // Helicity +1 -1
    TF1 *H_plus = new TF1("Hplus_draw", fitFunctor, xmin, xmax, 3);
    H_plus->SetParameters(f->GetParameters());
    H_plus->SetParameter(0, +1); 
    H_plus->SetLineColor(kBlue);
    H_plus->SetLineStyle(2);
    H_plus->SetLineWidth(2);
    H_plus->Draw("SAME");
    TF1 *H_minus = new TF1("Hminus_draw", fitFunctor, xmin, xmax, 3);
    H_minus->SetParameters(f->GetParameters());
    H_minus->SetParameter(0, -1); 
    H_minus->SetLineColor(kRed);
    H_minus->SetLineStyle(2);
    H_minus->SetLineWidth(2);
    H_minus->Draw("SAME");
    
    // legend
    TLegend* leg = new TLegend(0.15,0.15,0.35,0.35);
    leg->AddEntry(h, "signal", "l");
    leg->AddEntry(f, "global fit", "l");
    leg->AddEntry(H_plus, "H = +1", "l");
    leg->AddEntry(H_minus, "H = -1", "l");
    leg->Draw();

    gPad->Update();


	// file save
	const char* outdir =
        "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/MC/nofilters/";
        
	if (outname) {
		TString pdf = TString(outdir) + TString(outname) + "_analytic_fit.pdf";
		c->SaveAs(pdf);
		
	}

	
	delete h;
	delete f;
	delete c;
}

// FINAL ROOT VOID

void fitLep_analytic() {

    //gROOT->Reset();

    const char* infile =
        "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/histmaker/p8_ee_Ztautau_ecm91.root";

	TFile* f = TFile::Open(infile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: could not open file: " << infile << std::endl;
        return;
    }

    // Retrieve histograms 
    
    TH1D* h_el = (TH1D*)f->Get("el_sgn");
	TH1D* h_mu = (TH1D*)f->Get("mu_sgn");

	TH1D* h_MCel = (TH1D*)f->Get("mc_el_sgn");
	TH1D* h_MCmu = (TH1D*)f->Get("mc_mu_sgn");

	TH1D* h_el_free = (TH1D*)f->Get("el_free");
	TH1D* h_mu_free = (TH1D*)f->Get("mu_free");

    if (!h_el) std::cerr << "[WARN] Histogram el_sgn not found" << std::endl;
    if (!h_mu) std::cerr << "[WARN] Histogram mu_sgn not found" << std::endl;

    // Fit for electrons
    if (h_el) {
        TString el = "electron";
        TString title = "Reco events + Reco energy : e channel;x_{e};Events";
        Fit_1prong(h_el, el, title);
    }
	if (h_MCel) {
        TString el = "MCelectron";
        TString title = "Reco events + MC energy : e channel;x_{e};Events";
        Fit_1prong(h_MCel, el, title);
    }
	if (h_el_free) {
        TString el = "FREEelectron";
        TString title = "Reco NO invariant mass limit : e channel;x_{e};Events";
        Fit_1prong(h_el_free, el, title);
    }

    // Fit for muons
    if (h_mu) {
        TString mu = "muon";
        TString title = "Reco events + Reco energy : #mu channel;x_{#mu};Events";
        Fit_1prong(h_mu, mu, title);
    }
	if (h_MCmu) {
        TString mu = "MCmuon";
        TString title = "Reco events + MC energy : #mu channel;x_{#mu};Events";
        Fit_1prong(h_MCmu, mu, title);
    }
	if (h_mu_free) {
        TString mu  = "FREEmuon";
        TString title = "Reco NO invariant mass limit : #mu channel;x_{#mu};Events";
        Fit_1prong(h_mu_free, mu, title);
    }

    f->Close();
}


