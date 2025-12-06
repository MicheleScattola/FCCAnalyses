#include <iostream>
//================================================
double W_lep(double *x, double *par){
	double W =  par[1]*par[2]*1./3.*( (5-9*x[0]*x[0]+4*x[0]*x[0]*x[0]) + par[0] * (1-9*x[0]*x[0]+8*x[0]*x[0]*x[0]) );
	return W;
}

double aux_norm(double x, double P){
	return (5*3 - 3*x*x*x + x*x*x*x) + P * (x - 3*x*x*x + 2*x*x*x*x);
}
	

void Fit_1prong(TH1 *hist, const char *outname, const char *title) {
		
	if (!hist) {
		std::cerr << "[ERROR] Fit_1prong: hist is null\n";
		return;
	}

	// cloning and normalizing histogram
	TH1 *h = (TH1*)hist->Clone(Form("%s_clone", hist->GetName()));
	h->Sumw2();
	
	int N = h->Integral();
	cout << "\n =========================\n" << endl;
	cout << "[INFO] Fitting " << hist->GetName() << " with " << N << " entries." << endl;
	
	double bin_width = h->GetBinWidth(1);
	
	// TF1 with 1 parameter P
	//double xmin = h->GetXaxis()->GetXmin();
	//double xmax = h->GetXaxis()->GetXmax();
	double xmin = 0.05;
	double xmax = 1.0;
	TF1 *f = new TF1(Form("fit_%s", hist->GetName()),W_lep, xmin, xmax, 3);
	cout << "Fitting between " << xmin << " and " << xmax << endl;     
	// initialize values
	f->SetParName(0, "P_{#tau}");
	f->SetParameter(0, -0.1); 
	f->SetParameter(1, N);
	f->FixParameter(2, bin_width);
	
	double Perr = 0.0;
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1111);
	
	
	// plot
	TCanvas *c = new TCanvas(Form("canvas_%s", hist->GetName()), "P fit", 900,600);
	c->cd();
	
	// Fit:
	
	h->SetTitle(title);
	h->SetMarkerStyle(20);
	h->SetMarkerSize(0.7);
	h->SetLineColor(kBlack);
	h->SetLineWidth(1);
	      
	f->SetLineColor(kMagenta);
	f->SetLineWidth(2);
	
	h->SetStats(1);
	h->Draw("HIST E");
	h->Fit(f, "RQ SAME");
	f->Draw("SAME");
	
	// results
	double P = f->GetParameter(0);
	Perr = f->GetParError(0);

	std::cout << "[INFO] Fit obtained for " << hist->GetName() << ":" << endl;
	std::cout << "****** Polarization P = " << P << " +/- " << Perr << endl;
	cout << "\n =========================\n" << endl;
	
	// Additional drawings
    
    // Helicity +1 -1
    TF1 *H_plus = new TF1("Hplus_draw", W_lep, xmin, xmax, 3);
    H_plus->SetParameters(f->GetParameters());
    H_plus->SetParameter(0, +1); 
    H_plus->SetLineColor(kBlue);
    H_plus->SetLineStyle(2);
    H_plus->Draw("SAME");
    TF1 *H_minus = new TF1("Hminus_draw", W_lep, xmin, xmax, 3);
    H_minus->SetParameters(f->GetParameters());
    H_minus->SetParameter(0, -1); 
    H_minus->SetLineColor(kRed);
    H_minus->SetLineStyle(2);
    H_minus->Draw("SAME");
    
    // legend
    TLegend* leg = new TLegend(0.15,0.15,0.3,0.3);
    leg->AddEntry(h, "signal", "l");
    leg->AddEntry(f, "global fit", "l");
    leg->AddEntry(H_plus, "H = +1", "l");
    leg->AddEntry(H_minus, "H = -1", "l");
    leg->Draw();

    gPad->Update();


	// file save
	const char* outdir =
        "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/";
        
	if (outname) {
		TString png = TString(outdir) + TString(outname) + "_fit.png";
		TString pdf = TString(outdir) + TString(outname) + "_fit.pdf";
		c->SaveAs(png);
		c->SaveAs(pdf);
		
	}

	
	delete h;
	delete f;
	delete c;
}

// FINAL ROOT VOID

void fitLep() {

    //gROOT->Reset();

    const char* infile =
        "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/histmaker/p8_ee_Ztautau_ecm91.root";

    TFile* f = TFile::Open(infile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open file: " << infile << std::endl;
        return;
    }

    // histos
    TH1F* h_el = (TH1F*)f->Get("el_sgn");
    TH1F* h_mu = (TH1F*)f->Get("mu_sgn");

    if (!h_el) std::cerr << "[WARN] Histogram el_sgn not found" << std::endl;
    if (!h_mu) std::cerr << "[WARN] Histogram mu_sgn not found" << std::endl;

    // Fit for electrons
    if (h_el) {
        TString el = "electron";
        TString title = "#tau #rightarrow e #nu_{e} #nu_{#tau} channel;x_{e};Events";
        Fit_1prong(h_el, el, title);
    }

    // Fit for muons
    if (h_mu) {
        TString mu = "muon";
        TString title = "#tau #rightarrow #mu #nu_{#mu} #nu_{#tau} channel;x_{#mu};Events";
        Fit_1prong(h_mu, mu, title);
    }

    f->Close();
}


