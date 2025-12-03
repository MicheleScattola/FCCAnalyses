{
    gROOT->Reset();

    // Input and output paths (edit if needed)
    const char* infile = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/histmaker/p8_ee_Ztautau_ecm91.root";
    const char* outdir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/plots/confusion";

    TFile* f = TFile::Open(infile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: could not open file: " << infile << std::endl;
        return;
    }

    // Retrieve histograms 
    TH1F* rp_phe = (TH1F*)f->Get("RhoAsPi_phe");
    TH1F* rot_phe = (TH1F*)f->Get("RhoAsOt_phe");
    TH1F* otr_phe = (TH1F*)f->Get("OtAsRho_phe");
    
    TH1F* rp_pht = (TH1F*)f->Get("RhoAsPi_pht");
    TH1F* rot_pht = (TH1F*)f->Get("RhoAsOt_pht");
    TH1F* otr_pht = (TH1F*)f->Get("OtAsRho_pht");
    
    TH1F* rp_e = (TH1F*)f->Get("RhoAsPi_e");
    TH1F* ep_e = (TH1F*)f->Get("ElAsPi_e");
    TH1F* mp_e = (TH1F*)f->Get("MuAsPi_e");
    TH1F* pp_e = (TH1F*)f->Get("PiAsPi_e");
    TH1F* op_e = (TH1F*)f->Get("OtAsPi_e");
    
    TH1F* rp_t = (TH1F*)f->Get("RhoAsPi_t");
    TH1F* ep_t = (TH1F*)f->Get("ElAsPi_t");
    TH1F* mp_t = (TH1F*)f->Get("MuAsPi_t");
    TH1F* pp_t = (TH1F*)f->Get("PiAsPi_t");
    TH1F* op_t = (TH1F*)f->Get("OtAsPi_t");
    
    
    // Check histograms existence
    //if (!rp)   std::cerr << "ERROR: histogram 'rhoaspi' not found\n";
    //if (!ro)   std::cerr << "ERROR: histogram 'rhoasother' not found\n";
    //if (!otr)   std::cerr << "ERROR: histogram 'otherasrho' not found\n";


    // Crrpte output directory if missing (system call)
    gSystem->Exec(Form("mkdir -p %s", outdir));

    // -------------------------
    TCanvas* C = new TCanvas("C", "confusion Ph contributions", 1400, 500);
    //C->SetTitle("Miss-ID events: Photons contribution");
    C->Divide(2,1);
    C->SetGrid();
    C->cd(1);
    
    otr_phe->SetTitle("Reco #gamma energy spectrum;Energy_{#gamma} [GeV];Entries");
    otr_phe->SetLineColor(kOrange+7);
    otr_phe->SetLineWidth(2);
    otr_phe->Draw("HIST");
    
    rot_phe->SetLineColor(kGreen);
    rot_phe->SetLineWidth(2);
    rot_phe->Draw("HIST SAME");
    
    rp_phe->SetLineColor(kBlue);
    rp_phe->SetLineWidth(2);
    rp_phe->Draw("HIST SAME");
    
    // Legend
    TLegend* leg = new TLegend(0.60,0.60,0.88,0.88);
    leg->AddEntry(rp_phe, "TrueMC #rho as #pi", "l");
    leg->AddEntry(rot_phe, "TrueMC #rho as other", "l");
    leg->AddEntry(otr_phe, "TrueMC other as #rho", "l");
    leg->Draw();
    
    C->cd(2);
    
    otr_pht->SetTitle("Reco #gamma |cos#theta| ;cos#theta_{#gamma};Entries");
    otr_pht->SetLineColor(kOrange+7);
    otr_pht->SetLineWidth(2);
    otr_pht->Draw("HIST");
    
    rot_pht->SetLineColor(kGreen);
    rot_pht->SetLineWidth(2);
    rot_pht->Draw("HIST SAME");
    
    rp_pht->SetLineColor(kBlue);
    rp_pht->SetLineWidth(2);
    rp_pht->Draw("HIST SAME");
    
    // Legend
    TLegend* legg = new TLegend(0.40,0.60,0.70,0.88);
    legg->AddEntry(rp_pht, "TrueMC #rho as #pi", "l");
    legg->AddEntry(rot_pht, "TrueMC #rho as other", "l");
    legg->AddEntry(otr_pht, "TrueMC other as #rho", "l");
    legg->Draw();

    C->Update();
    C->SaveAs(Form("%s/miss_id_ph.png", outdir));
    C->SaveAs(Form("%s/miss_id_ph.pdf", outdir));
    
    // -------------------------
    // Single canvas with all plots
    TCanvas* C1 = new TCanvas("C1", "confusion 1-prong pi", 1400, 500);
    //C1->SetTitle("Miss-ID events: Reco #pi contribution");
    C1->Divide(2,1);
    C1->SetGrid();
    C1->cd(1);
    
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
    
    op_e->SetLineColor(kGray+3);
    op_e->SetLineWidth(2);
    op_e->Draw("HIST SAME");
    
    // Legend
    TLegend* leg1 = new TLegend(0.60,0.60,0.88,0.88);
    leg1->AddEntry(pp_e, "TrueMC #pi as #pi", "l");
    leg1->AddEntry(rp_e, "TrueMC #rho as #pi", "l");
    leg1->AddEntry(mp_e, "TrueMC #mu as #pi", "l");
    leg1->AddEntry(ep_e, "TrueMC e^{#pm} as #pi", "l");
    leg1->AddEntry(op_e, "TrueMC other as #pi", "l");
    leg1->Draw();
    
    C1->cd(2);
    
    rp_t->SetTitle("Reco #pi |cos#theta| ;cos#theta_{#pi};Entries");
    rp_t->SetLineColor(kRed);
    rp_t->SetLineWidth(2);
    rp_t->Draw("HIST");
    
    mp_t->SetLineColor(kMagenta);
    mp_t->SetLineWidth(2);
    mp_t->Draw("HIST SAME");
    
    ep_t->SetLineColor(kBlue);
    ep_t->SetLineWidth(2);
    ep_t->Draw("HIST SAME");
    
    pp_t->SetLineColor(kGreen+2);
    pp_t->SetLineWidth(2);
    pp_t->Draw("HIST SAME");
    
    op_t->SetLineColor(kGray+3);
    op_t->SetLineWidth(2);
    op_t->Draw("HIST SAME");
    
    // Legend
    TLegend* legg1 = new TLegend(0.40,0.60,0.70,0.88);
    legg1->AddEntry(pp_t, "TrueMC #pi as #pi", "l");
    legg1->AddEntry(rp_t, "TrueMC #rho as #pi", "l");
    legg1->AddEntry(mp_t, "TrueMC #mu as #pi", "l");
    legg1->AddEntry(ep_t, "TrueMC e^{#pm} as #pi", "l");
    legg1->AddEntry(op_t, "TrueMC other as #pi", "l");
    legg1->Draw();
    
    C1->Update();
    C1->SaveAs(Form("%s/miss_id_pi1prong.png", outdir));
    C1->SaveAs(Form("%s/miss_id_pi1prong.pdf", outdir));
    
    
    
    // Clean up 
    // f->Close();
    // delete f;
}


