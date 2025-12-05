{
    gROOT->Reset();

    // Input and output paths (edit if needed)
    const char* infile = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau//histmaker/p8_ee_Ztautau_ecm91.root";
    const char* outdir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/plots/confusion";

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
    
    /*
    TH1F* rp_t = (TH1F*)f->Get("RhoAsPi_t");
    TH1F* ep_t = (TH1F*)f->Get("ElAsPi_t");
    TH1F* mp_t = (TH1F*)f->Get("MuAsPi_t");
    TH1F* pp_t = (TH1F*)f->Get("PiAsPi_t");
    TH1F* op_t = (TH1F*)f->Get("OtAsPi_t");
    TH1F* ap_t = (TH1F*)f->Get("A1AsPi_t");
    */
    
    
    // Check histograms existence

    // Create output directory if missing (system call)
    gSystem->Exec(Form("mkdir -p %s", outdir));
    
    // -------------------------
    TCanvas* C1 = new TCanvas("C1", "confusion 1-prong pi", 800, 600);
    //C1->SetTitle("Miss-ID events: Reco #pi contribution");
    C1->SetGrid();
    C1->cd();
    TGaxis::SetMaxDigits(3);
    
    mp_e->SetTitle("Reco #pi energy spectrum;Energy_{#pi} [GeV];Entries");
    mp_e->SetLineColor(kMagenta);
    mp_e->SetLineWidth(2);
    mp_e->Draw("HIST");
    
    ep_e->SetLineColor(kRed);
    ep_e->SetLineWidth(2);
    ep_e->Draw("HIST SAME");
    
    pp_e->SetLineColor(kGreen+2);
    pp_e->SetLineWidth(2);
    pp_e->Draw("HIST SAME");
    
    rp_e->SetLineColor(kBlue);
    rp_e->SetLineWidth(2);
    rp_e->Draw("HIST SAME");

    ap_e->SetLineColor(kAzure+10);
    ap_e->SetLineWidth(2);
    ap_e->Draw("HIST SAME");
    
    op_e->SetLineColor(kGray+3);
    op_e->SetLineWidth(2);
    op_e->SetFillColor(kGray+3); 
    op_e->Draw("HIST SAME");
    
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
    

    /*
    TCanvas* C2 = new TCanvas("C2", "confusion 1-prong pi - costheta", 800, 600);
    //C2->SetTitle("Miss-ID events: Reco #pi contribution");
    C2->SetGrid();
    C2->cd();
    TGaxis::SetMaxDigits(3);
    
    rp_t->SetTitle("Reco #pi |cos#theta| ;cos#theta_{#pi};Entries");
    rp_t->SetLineColor(kRed);
    rp_t->SetLineWidth(2);
    rp_t->Draw("HIST");
    
    mp_t->SetLineColor(kMagenta);
    mp_t->Add(ep_t);
    mp_t->SetLineWidth(2);
    mp_t->Draw("HIST SAME");
    
    //ep_t->SetLineColor(kRed);
    //ep_t->SetLineWidth(2);
    //ep_t->Draw("HIST SAME");
    
    pp_t->SetLineColor(kGreen+2);
    pp_t->SetLineWidth(2);
    pp_t->Draw("HIST SAME");

    ap_t->SetLineColor(kAzure+10);
    ap_t->SetLineWidth(2);
    ap_t->Draw("HIST SAME");
    
    op_t->SetLineColor(kGray+3);
    op_t->SetLineWidth(2);
    op_t->SetFillColor(kGray+3);
    op_t->Draw("HIST SAME");
    
    // Legend
    TLegend* legg1 = new TLegend(0.40,0.60,0.70,0.88);
    legg1->AddEntry(mp_t, "TrueMC lepton as #pi", "l");
    //legg1->AddEntry(ep_t, "TrueMC e^{#pm} as #pi", "l");
    legg1->AddEntry(pp_t, "TrueMC #pi as #pi", "l");
    legg1->AddEntry(rp_t, "TrueMC #rho as #pi", "l");
    legg1->AddEntry(ap_t, "TrueMC #rho as #pi", "l");
    legg1->AddEntry(op_t, "TrueMC other as #pi", "l");
    legg1->Draw();
    
    C2->Update();
    C2->SaveAs(Form("%s/missid_costheta.png", outdir));
    C2->SaveAs(Form("%s/missid_costheta.pdf", outdir));
    */
    
    
    // Clean up 
    // f->Close();
    // delete f;
}


