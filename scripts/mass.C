void mass(){

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
    
    TH1F* aMC = (TH1F*)f->Get("MC_a1_m");
    TH1F* rMC = (TH1F*)f->Get("MC_rho_m");
    TH1F* aReco = (TH1F*)f->Get("reco_a1_m");
    TH1F* rReco = (TH1F*)f->Get("reco_rho_m");
    
    gSystem->Exec(Form("mkdir -p %s", outdir));
    
    // -------------------------
    TCanvas* C1 = new TCanvas("C1", "MC masses", 800, 600);
    C1->SetGrid();
    C1->cd();
    TGaxis::SetMaxDigits(3);
    
    aMC->SetTitle("MC invariant masses;Mass [GeV];Entries");
    aMC->SetLineColor(kBlue);
    aMC->SetFillColorAlpha(kBlue,0.3);
    aMC->SetLineWidth(2);
    aMC->Draw("HIST");
    
    rMC->SetLineColor(kRed);
    rMC->SetFillColorAlpha(kRed,0.3);
    rMC->SetLineWidth(2);
    rMC->Draw("HIST SAME");
    
    // Legend
    TLegend* leg1 = new TLegend(0.65,0.70,0.88,0.88);
    leg1->AddEntry(aMC, "TrueMC a_{1}", "l");
    leg1->AddEntry(rMC, "TrueMC #rho", "l");
    leg1->Draw();

    C1->Update();
    C1->SaveAs(Form("%s/MC_masses.pdf", outdir));

    // ============================
    // now same but with reco masses
    TCanvas* C2 = new TCanvas("C2", "Reco masses", 800, 600);
    C2->SetGrid();
    C2->cd();
    TGaxis::SetMaxDigits(3);

    aReco->SetTitle("Reconstructed invariant masses;Mass [GeV];Entries");
    aReco->SetLineColor(kBlue);
    aReco->SetFillColorAlpha(kBlue,0.3);
    aReco->SetLineWidth(2);
    aReco->Draw("HIST");    

    rReco->SetLineColor(kRed);
    rReco->SetFillColorAlpha(kRed,0.3);
    rReco->SetLineWidth(2);
    rReco->Draw("HIST SAME");
    
    // Legend
    TLegend* leg2 = new TLegend(0.65,0.70,0.88,0.88);
    leg2->AddEntry(aMC, "Reco a_{1}", "l");
    leg2->AddEntry(rMC, "Reco #rho", "l");
    leg2->Draw();

    C1->Update();
    C1->SaveAs(Form("%s/Reco_masses.pdf", outdir));

    // ============================
    // now pulls
    TCanvas* C3 = new TCanvas("C3", "Mass pulls", 800, 600);
    C3->Divide(2,1);
    C3->cd(1);
    C3->cd(1)->SetGrid();
    TGaxis::SetMaxDigits(3);
    TH1F* aPull = (TH1F*)aMC->Clone("aPull");
    aPull->Add(aReco, -1);
    aPull->SetTitle("a_{1} mass pull;Mass_{MC} - Mass_{Reco} [GeV];Entries");
    aPull->SetLineColor(kBlue);
    aPull->SetFillColorAlpha(kBlue,0.3);
    aPull->SetLineWidth(2);
    aPull->Draw("HIST");

    C3->cd(2);
    C3->cd(2)->SetGrid();
    TGaxis::SetMaxDigits(3);
    TH1F* rPull = (TH1F*)rMC->Clone("rPull");
    rPull->Add(rReco, -1);
    rPull->SetTitle("#rho mass pull;Mass_{MC} - Mass_{Reco} [GeV];Entries");
    rPull->SetLineColor(kRed);
    rPull->SetFillColorAlpha(kRed,0.3);
    rPull->SetLineWidth(2);
    rPull->Draw("HIST");

    C3->Update();
    C3->SaveAs(Form("%s/mass_pulls.pdf", outdir));

    // Clean up 
    // f->Close();
    // delete f;
}


