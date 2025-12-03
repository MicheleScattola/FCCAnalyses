#include <cmath>
//================================================
double BKG(double *x, double *par){
	// gaussian + linear noise
	return (par[0]/(par[2]*sqrt(2*M_PI)))*exp(-(x[0]-par[1])*(x[0]-par[1])/(2*par[2]*par[2])) + par[3]*x[0] + par[4];
}

void Fit_bkg(TH1 *hist) {
		
	if (!hist) {
		std::cerr << "[ERROR] Fit_Pi_Background: hist is null\n";
		return;
	}

	// cloning and normalizing histogram
	TH1 *h = (TH1*)hist->Clone(Form("%s_clone", hist->GetName()));
	double integral = h->Integral("width");
	int N = h->Integral();
	if (integral <= 0) {
		std::cerr << "[ERROR] Fit_Pi_Background: hist integral (width) <= 0\n";
		delete h;
		return;
	} else {
		cout << "\n =========================\n" << endl;
		cout << "[INFO] Fitting " << hist->GetName() << " with " << N << " entries." << endl;}
	h->Scale(1.0 / integral);

	// TF1 with 1 parameter P
	//double xmin = h->GetXaxis()->GetXmin();
	//double xmax = h->GetXaxis()->GetXmax();
	double xmin = 0.05;
	double xmax = 0.95;
	TF1 *f = new TF1("bkg_shape",BKG, xmin, xmax, 5);
	cout << "Fitting between " << xmin << " and " << xmax << endl;     
	// initialize values
	f->SetParNames("N","x_{0}","#sigma","a","b");
	f->SetParameter(0, 0.2);  // N (gaussia area)
    f->SetParameter(1, 0.2);  // mean
    f->SetParameter(2, 0.1);  // sigma 
    f->SetParameter(3, -0.7); // slope
    f->SetParameter(4, 1.2);  // offset
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0111);
	
	
	// plot
	TCanvas *c = new TCanvas(Form("canvas_%s", hist->GetName()), "bkg fit", 800,600);
	c->cd();
	
	// Fit:
	
	h->SetTitle("Background fit (gaus + linear);x;1/#Gamma d#Gamma/dx");
	h->SetMarkerStyle(20);
	h->SetMarkerSize(0.7);
	h->SetLineColor(kBlack);
	h->SetLineWidth(1);
	      
	f->SetLineColor(kRed);
	f->SetLineWidth(2);
	
	h->SetStats(1);
	h->Draw("HIST E");
	h->Fit(f, "RQ SAME");
	f->Draw("SAME");
	
	// results

	cout << "[INFO] Fit obtained for " << hist->GetName() << ":" << endl;
	cout << "****** Background estimated with gaussian + linear noise:" << endl;
	cout << "****** Gaussian: N = " << f->GetParameter(0) << " | x_{0} = " << f->GetParameter(1) << " | sigma = " << f->GetParameter(2) << endl;
	cout << "****** Linear: " << f->GetParameter(3) << " x + " << f->GetParameter(4) << endl;
	cout << "\n =========================\n" << endl;

    gPad->Update();


	// file save
	const char* outdir =
        "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/plots/optimal/";
        
	
	TString png = TString(outdir) + "piBKG_gauslin.png";
	TString pdf = TString(outdir) + "piBKG_gauslin.pdf";
	c->SaveAs(png);
	c->SaveAs(pdf);
		

	
	//delete h;
	//delete c;
}

// FINAL ROOT VOID

void bkg_gauslin() {

    //gROOT->Reset();

    const char* infile =
        "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/histmaker/optimal/p8_ee_Ztautau_ecm91.root";

    TFile* f = TFile::Open(infile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open file: " << infile << std::endl;
        return;
    }

    // histos
    TH1F* h = (TH1F*)f->Get("pi_bkg");

    if (!h) cerr << "[WARN] Histogram el_sgn not found" << std::endl;

    Fit_bkg(h);

    //f->Close();
}


