#include "ROOT/TMath.cxx"

// invariant mass BreitWigner distribution
float BreitWigner(double *x, double *par){
	return TMath::BreitWigner (x[0],par[0],par[1]);
}

void mass() {
    
    ROOT::EnableImplicitMT(); 

    gROOT->Reset();

    const char* infile = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/treemaker/p8_ee_Ztautau_ecm91.root";
    const char* outdir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/";

    gSystem->Exec(Form("mkdir -p %s", outdir));

    ROOT::RDataFrame df("events", infile);

    int nbins_m = 60; double min_m = 0; double max_m = 2;
    int nbins_p = 50; double min_p = -0.6; double max_p = 0.6;

    // MC Masses
    auto h_aMC = df.Histo1D({"MC_a1_m", "MC a_{1} mass", nbins_m, min_m, max_m}, "MC_a1_m");
    auto h_rMC = df.Histo1D({"MC_rho_m", "MC #rho mass", nbins_m, min_m, max_m}, "MC_rho_m");
    
    // Reco Masses
    auto h_aReco = df.Histo1D({"reco_a1_m", "Reconstructed a_{1} mass", nbins_m, min_m, max_m}, "reco_a1_m");
    auto h_rReco = df.Histo1D({"reco_rho_m", "Reconstructed #rho mass", nbins_m, min_m, max_m}, "reco_rho_m");
    
    // Pulls
    auto h_aPull = df.Histo1D({"a1_pull", "a_{1} mass pull", nbins_p, min_p, max_p}, "a1_pull");
    auto h_rPull = df.Histo1D({"rho_pull", "#rho mass pull", nbins_p, min_p, max_p}, "rho_pull");

    // Clone histograms to detach from RDataFrame
    TH1D* aMC = (TH1D*)h_aMC->Clone();
    TH1D* rMC = (TH1D*)h_rMC->Clone();
    TH1D* aReco = (TH1D*)h_aReco->Clone();
    TH1D* rReco = (TH1D*)h_rReco->Clone();
    TH1D* aPull = (TH1D*)h_aPull->Clone();
    TH1D* rPull = (TH1D*)h_rPull->Clone();

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1111);
    // A1
    TCanvas* C1 = new TCanvas("C1", "a1 Comparison", 800, 600);
    C1->SetGrid();
    C1->cd();
    TGaxis::SetMaxDigits(3);

    aMC->SetTitle("a_{1} Mass: MC vs Reco;Mass [GeV];Entries");
    aMC->SetLineColor(kBlue + 1);
    aMC->SetFillColorAlpha(kBlue + 1, 0.3); 
    aMC->SetLineWidth(2);
    aMC->Draw("HIST E");
    
    aReco->SetLineColor(kRed);
    aReco->SetLineWidth(2);
    // aReco->SetFillColorAlpha(kRed, 0.0);
    // convolution fit
    TF1* f = new TF1("f", BreitWigner, 0, 2, 2);
    f->SetParameters(1.26, 0.4); // initial values
    f->SetParNames("mass", "#Gamma");
    aReco->Draw("HIST SAMES E");
    aReco->Fit(f,"RQ SAME");
    f->Draw("SAME");
    C1->Update();
    
    TLegend* leg1 = new TLegend(0.2, 0.70, 0.4, 0.88);
    leg1->AddEntry(aMC, "MC a_{1}", "f");   
    leg1->AddEntry(aReco, "Reco a_{1}", "l");
    leg1->Draw();

    C1->SaveAs(Form("%s/a1_mass.pdf", outdir));

    // RHO
    TCanvas* C2 = new TCanvas("C2", "rho Comparison", 800, 600);
    C2->SetGrid();
    C2->cd();
    TGaxis::SetMaxDigits(3);

    rMC->SetTitle("#rho Mass: MC vs Reco;Mass [GeV];Entries");
    rMC->SetLineColor(kBlue + 1);
    rMC->SetFillColorAlpha(kBlue + 1, 0.3); 
    rMC->SetLineWidth(2);
    rMC->Draw("HIST E");
    
    rReco->SetLineColor(kRed);
    rReco->SetLineWidth(2);
    rReco->Fit("gaus","Q"); 
    rReco->Draw("HIST SAME S E");   
    C2->Update();
    
    TLegend* leg2 = new TLegend(0.65, 0.30, 0.88, 0.48);
    leg2->AddEntry(rMC, "MC #rho", "f");
    leg2->AddEntry(rReco, "Reco #rho", "l");
    leg2->Draw();

    C2->SaveAs(Form("%s/rho_mass.pdf", outdir));

    // PULLS
    TCanvas* C3 = new TCanvas("C3", "Mass pulls", 900, 600);
    C3->Divide(2, 1);
    
    C3->cd(1);
    C3->cd(1)->SetGrid();
    TGaxis::SetMaxDigits(3);
    aPull->SetTitle("a_{1} mass pull;Mass_{MC} - Mass_{Reco} [GeV];Entries");
    aPull->SetLineColor(kBlue+1);
    aPull->SetFillColorAlpha(kBlue, 0.3);
    aPull->SetLineWidth(2);
    aPull->Draw("HIST");

    C3->cd(2);
    C3->cd(2)->SetGrid();
    TGaxis::SetMaxDigits(3);
    rPull->SetTitle("#rho mass pull;Mass_{MC} - Mass_{Reco} [GeV];Entries");
    rPull->SetLineColor(kRed);
    rPull->SetFillColorAlpha(kRed, 0.3);
    rPull->SetLineWidth(2);
    rPull->Draw("HIST");

    C3->SaveAs(Form("%s/mass_pulls.pdf", outdir));
    
    // ==========================================
    // OVERFLOW COUNTS
    // ==========================================
    
    // Find bin for 2 GeV
    int bin_a_start = aReco->FindBin(2);
    int bin_r_start = rReco->FindBin(2);
    
    // GetNbinsX() + 1 is the OVERFLOW bin
    double counts_a = aReco->Integral(bin_a_start, aReco->GetNbinsX() + 1);
    double counts_r = rReco->Integral(bin_r_start, rReco->GetNbinsX() + 1);

    std::cout << "\n=======================================" << std::endl;
    std::cout << "COUNTS FOR RECO MASS > 2 GeV" << std::endl;
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "Reco a1 counts : " << counts_a << std::endl;
    std::cout << "Reco rho counts: " << counts_r << std::endl;
    std::cout << "TOTAL          : " << counts_a + counts_r << std::endl;
    std::cout << "=======================================\n" << std::endl;
}
