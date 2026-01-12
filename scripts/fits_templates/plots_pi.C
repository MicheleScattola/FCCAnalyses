void plots_pi(){
    gROOT->Reset();

    // Input and output paths (edit if needed)
    const char* infile = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/histmaker/p8_ee_Ztautau_ecm91.root";
    const char* outdir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/plots/";

    TFile* f = TFile::Open(infile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: could not open file: " << infile << std::endl;
        return;
    }

    // Retrieve histograms 
    
    TH1F* rp_e = (TH1F*)f->Get("RhoAsPi_e");
    TH1F* ep_e = (TH1F*)f->Get("ElAsPi_e");
    TH1F* mp_e = (TH1F*)f->Get("MuAsPi_e");
    TH1F* pp_e = (TH1F*)f->Get("PiAsPi_e");
    TH1F* op_e = (TH1F*)f->Get("OtAsPi_e");
    TH1F* ap_e = (TH1F*)f->Get("A1AsPi_e");
    
    gSystem->Exec(Form("mkdir -p %s", outdir));
    
    // -------------------------
    TCanvas* C1 = new TCanvas("C1", "confusion 1-prong pi", 800, 600);
    //C1->SetTitle("Miss-ID events: Reco #pi contribution");
    C1->SetGrid();
    C1->cd();
    TGaxis::SetMaxDigits(3);
    
    mp_e->SetTitle("Reco #pi energy spectrum;Energy_{#pi} [GeV];Entries");
    mp_e->SetLineColor(kMagenta);
    mp_e->SetFillColorAlpha(kMagenta,0.3);
    mp_e->SetLineWidth(2);
    mp_e->Draw("HIST");
    
    ep_e->SetLineColor(kRed);
    ep_e->SetFillColorAlpha(kRed,0.3);
    ep_e->SetLineWidth(2);
    ep_e->Draw("HIST SAME");
    
    pp_e->SetLineColor(kGreen+2);
    pp_e->SetLineWidth(2);
    pp_e->SetFillColorAlpha(kGreen+2,0.3);
    pp_e->Draw("HIST SAME");
    
    rp_e->SetLineColor(kBlue);
    rp_e->SetLineWidth(2);
    rp_e->SetFillColorAlpha(kBlue,0.3);
    rp_e->Draw("HIST SAME");
    
    op_e->SetLineColor(kGray+3);
    op_e->SetLineWidth(2);
    op_e->SetFillColor(kGray+3); 
    op_e->Draw("HIST SAME");
    
    ap_e->SetLineColor(kAzure+10);
    ap_e->SetLineWidth(2);
    ap_e->SetFillColor(kAzure+10);
    ap_e->Draw("HIST SAME");
    
    // Legend
    TLegend* leg1 = new TLegend(0.60,0.60,0.88,0.88);
    leg1->AddEntry(mp_e, "TrueMC #mu as #pi", "l");
    leg1->AddEntry(ep_e, "TrueMC e^{#pm} as #pi", "l");
    leg1->AddEntry(pp_e, "TrueMC #pi as #pi", "l");
    leg1->AddEntry(rp_e, "TrueMC #rho as #pi", "l");
    leg1->AddEntry(ap_e, "TrueMC a_{1} as #pi", "l");
    leg1->AddEntry(op_e, "TrueMC other as #pi", "l");
    leg1->Draw();

    C1->Update();
    C1->SaveAs(Form("%s/missid_e.png", outdir));
    C1->SaveAs(Form("%s/missid_e.pdf", outdir));
    

    
    // Clean up 
    // f->Close();
    // delete f;
}


