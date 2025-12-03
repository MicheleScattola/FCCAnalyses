{
    gROOT->Reset();

    // Input and output paths (edit if needed)
    const char* infile = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/histmaker/p8_ee_Ztautau_ecm91.root";
    const char* outdir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/plots/rho";

    TFile* f = TFile::Open(infile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: could not open file: " << infile << std::endl;
        return;
    }

    // Retrieve histograms 
    TH1F* gg = (TH1F*)f->Get("gammagamma_e");
    TH1F* ph = (TH1F*)f->Get("photons_e");
    TH1F* MCgg = (TH1F*)f->Get("MC_gammagamma_e");
    
    TH1F* p1g1t = (TH1F*)f->Get("pi1gamma1_true");
    TH1F* p1g1f = (TH1F*)f->Get("pi1gamma1_false");
    TH1F* p1g2t = (TH1F*)f->Get("pi1gamma2_true");
    TH1F* p1g2f = (TH1F*)f->Get("pi1gamma2_false");
    
    TH1F* w1p = (TH1F*)f->Get("wrongrho_1ph");
    TH1F* w2p = (TH1F*)f->Get("wrongrho_2ph");
    
    TH1F* sum2f = (TH1F*)f->Get("sum_wrong2ph");
    TH1F* sum2t = (TH1F*)f->Get("sum_true2ph");
    
    TH1F* tau = (TH1F*)f->Get("taus_e_collinear");
    
    TH1F* pdg1 = (TH1F*)f->Get("parentPDG1");
    TH1F* pdg2 = (TH1F*)f->Get("parentPDG2");


    // Create output directory if missing (system call)
    gSystem->Exec(Form("mkdir -p %s", outdir));

    // -------------------------
    // Single canvas with all plots
    TCanvas* C = new TCanvas("C", "MC vs Reco photons", 800, 500);
    C->SetGrid();
    C->cd();

    // Draw Reco2 gammagamma first
    if (gg) {
        gg->SetLineColor(kOrange+10);
        gg->SetLineWidth(2);
        gg->SetTitle("MC2Reco_energy vs MC_energy : #rho #rightarrow #gamma #gamma vs AllReco #gamma;Energy_{#gamma} [GeV];Entries");
        gg->Draw("HIST");
    } else {
        TLatex t; t.SetNDC(); t.DrawLatex(0.2,0.6,"Missing histogram: gammagamma_e");
    }

    // Overlay generic reconstructed photons
    if (ph) {
        ph->SetLineColor(kAzure);
        ph->SetLineWidth(2);
        if (gg) {
            // Ensure y-range accommodates both
            double ymax = std::max(gg->GetMaximum(), ph->GetMaximum());
            gg->SetMaximum(1.1*ymax);
        }
        ph->Draw("HIST SAME");
    } else {
        TLatex t; t.SetNDC(); t.DrawLatex(0.2,0.5,"Missing histogram: photons_e");
    }
    
    // Overlay MC true energy
    if (MCgg) {
        MCgg->SetLineColor(kViolet+1);
        MCgg->SetLineWidth(2);
        if (gg) {
            // Ensure y-range accommodates both
            double ymax = std::max(gg->GetMaximum(), ph->GetMaximum());
            gg->SetMaximum(1.1*ymax);
        }
        MCgg->Draw("HIST SAME");
    } else {
        TLatex t; t.SetNDC(); t.DrawLatex(0.2,0.5,"Missing histogram: MC_gammagamma_e");
    }
    
    
    // Legend
    TLegend* leg = new TLegend(0.60,0.75,0.88,0.88);
    if (ph) leg->AddEntry(ph, "All reconstructed #gamma", "l");
    if (gg) leg->AddEntry(gg, "#rho #rightarrow #gamma#gamma (MC2Reco)", "l");
    if (MCgg) leg->AddEntry(MCgg, "#rho #rightarrow #gamma#gamma (TrueMC)", "l");
    leg->Draw();

    C->Update();
    C->SaveAs(Form("%s/MC_vs_RP_gamma.png", outdir));
    C->SaveAs(Form("%s/MC_vs_RP_gamma.pdf", outdir));
    
    // -------------------------
    // Canvas 2: 
    TCanvas* C2 = new TCanvas("C2", "1 pi + n gamma", 800, 500);
    C2->SetGrid();
    C2->cd();
    
    if (p1g2f) {
        p1g2f->SetLineColor(kRed+1);
        p1g2f->SetLineWidth(2);
        p1g2f->SetTitle("Reconstructed 1 #pi + N #gamma : TrueMC vs ALL events;Energy_{#gamma} [GeV];Entries");
        p1g2f->Draw("HIST");
    }

    if (p1g2t) {
        p1g2t->SetLineColor(kMagenta);
        p1g2t->SetLineWidth(2);
        p1g2t->Draw("HIST SAME");
    }
    
    if (p1g1f) {
        p1g1f->SetLineColor(kBlue+1);
        p1g1f->SetLineWidth(2);
        p1g1f->Draw("HIST SAME");
    }

    if (p1g1t) {
        p1g1t->SetLineColor(kCyan);
        p1g1t->SetLineWidth(2);
        p1g1t->Draw("HIST SAME");
    }


    TLegend* leg2 = new TLegend(0.60,0.75,0.88,0.88);
    if (p1g1t) leg2->AddEntry(p1g1t, "Reco 1#pi + 1#gamma | True #rho evts", "l");
    if (p1g1f) leg2->AddEntry(p1g1f, "Reco 1#pi + 1#gamma | ALL events", "l");
    if (p1g2t) leg2->AddEntry(p1g2t, "Reco 1#pi + 2#gamma | True #rho evts", "l");
    if (p1g2f) leg2->AddEntry(p1g2f, "Reco 1#pi + 2#gamma | ALL events", "l");
    leg2->Draw();
    
    // zoom pad
	TPad *pad_zoom = new TPad("pad_zoom","zoom",0.40,0.2,0.8,0.7);
	pad_zoom->SetFillColor(0);
	pad_zoom->SetLineColor(1);
	pad_zoom->SetGrid();
	pad_zoom->Draw();
	pad_zoom->cd();

	// cloning
	TH1F *h1 = (TH1F*)p1g1t->Clone("p1g1t_zoom");
	TH1F *h2 = (TH1F*)p1g1f->Clone("p1g1f_zoom");
	TH1F *h3 = (TH1F*)p1g2t->Clone("p1g2t_zoom");
	TH1F *h4 = (TH1F*)p1g2f->Clone("p1g2f_zoom");

	// range
	h1->GetXaxis()->SetRangeUser(20.,35.);
	h1->GetYaxis()->SetRangeUser(0.,3500.);
	
	h1->SetTitle(";Energy_{#gamma} [GeV];");

	h1->Draw("HIST");
	h2->Draw("HIST SAME");
	h3->Draw("HIST SAME");
	h4->Draw("HIST SAME");
	
	C2->cd();

    C2->SaveAs(Form("%s/gamma_compare.png", outdir));
    C2->SaveAs(Form("%s/gamma_compare.pdf", outdir));
    
    // -------------------------
    // Canvas 3: 
    TCanvas* C3 = new TCanvas("C3", "NON-rho events 1 pi + n gamma", 800, 500);
    C3->SetGrid();
    C3->cd();
    
    if (w2p) {
        w2p->SetLineColor(kRed+1);
        w2p->SetLineWidth(2);
        w2p->SetTitle("Reconstructed 1 #pi + N #gamma : wrong non true #rho events;Energy_{#gamma} [GeV];Entries");
        w2p->Draw("HIST");
    }
    
    if (w1p) {
        w1p->SetLineColor(kBlue+1);
        w1p->SetLineWidth(2);
        w1p->Draw("HIST SAME");
    }
    
   

    TLegend* leg3 = new TLegend(0.60,0.75,0.88,0.88);
    if (w1p) leg3->AddEntry(w1p, "Reco 1#pi + 1#gamma | NOT a #rho", "l");
    if (w2p) leg3->AddEntry(w2p, "Reco 1#pi + 2#gamma | NOT a #rho", "l");
    //if (p1g1t) leg3->AddEntry(p1g1t, "Reco 1#pi + 1#gamma | TrueMC #rho", "l");
    leg3->Draw();

    C3->SaveAs(Form("%s/non_rho.png", outdir));
    C3->SaveAs(Form("%s/non_rho.pdf", outdir));
    
    // Canvas 4: 
    TCanvas* C4 = new TCanvas("C4", "test2", 800, 500);
    C4->SetGrid();
    C4->cd();
    
    sum2t->SetLineColor(kMagenta);
    sum2t->SetLineWidth(2);
    sum2t->SetTitle("Reco 1 #pi + 2 #gamma : SUM of #gamma energies;Energy_{#gamma} [GeV];counts");
    sum2t->Draw("HIST");
    
    sum2f->SetLineColor(kRed+1);
    sum2f->SetLineWidth(2);
    sum2f->Draw("HIST SAME");
    
    TLegend* leg5 = new TLegend(0.60,0.75,0.88,0.88);
	leg5->AddEntry(sum2f, "NON #rho events", "l");
	leg5->AddEntry(sum2t, "TrueMC #rho events", "l");
    leg5->Draw();
    
    C4->SaveAs(Form("%s/sum2ph.png", outdir));
    C4->SaveAs(Form("%s/sum2ph.pdf", outdir));
    
    // Canvas 5: 
    TCanvas* C5 = new TCanvas("C5", "ratio", 800, 500);
    C5->SetGrid();
    C5->cd();
    
	
	//Rapporto
	TH1F* ratio1 = (TH1F*)w1p->Clone("ratio1");
	ratio1->Divide(p1g1f);              // 1pi 1gamma: ratio of wrong non-rho in all Reco with 1pi1gamma
	TH1F* ratio2 = (TH1F*)w2p->Clone("ratio2");
	ratio2->Divide(p1g2f);              // 1pi 2gamma: ratio of wrong non-rho in all Reco with 1pi1gamma
	
	TH1F* ratio3 = (TH1F*)sum2f->Clone("ratio3");
	ratio3->Divide(sum2t);              // 1pi 2gamma: ratio of wrong nonpi-0 in SUM of the 2gammas energies
	
	ratio1->SetTitle("Ratio of wrong NON #rho events in ALL Reconstructed 1#pi + N#gamma;Energy_{#gamma} [GeV];Ratio");
	
	ratio2->SetLineColor(kOrange+10);
	ratio1->SetLineColor(kAzure);
	ratio3->SetLineColor(kGreen+2);
	
	ratio1->Draw("HIST");
	ratio2->Draw("HIST SAME");
	//ratio3->Draw("HIST SAME");
	
	TLegend* leg4 = new TLegend(0.15,0.65,0.4,0.88);
    leg4->AddEntry(ratio1, "Ratio in 1#pi + 1#gamma evts", "l");
    leg4->AddEntry(ratio2, "Ratio in 1#pi + 2#gamma evts", "l");
    //leg4->AddEntry(ratio3, "Ratio for SUM E_{1}+E_{2}", "l");
    leg4->Draw();
    
    C5->SaveAs(Form("%s/ratio.png", outdir));
    C5->SaveAs(Form("%s/ratio.pdf", outdir));
    
    // Canvas 6: 
    TCanvas* C6 = new TCanvas("C6", "taus", 800, 500);
    C6->SetGrid();
    C6->cd();
    
    tau->SetLineColor(kBlue);
    tau->SetLineWidth(2);
    tau->SetTitle("#tau reconstructed energies : collinear approximation;Energy_{#tau} [GeV];Entries");
    //gStyle->SetOptStat(1110);
	gStyle->SetOptFit(1111);
	//tau->SetStats(1);
	tau->Fit("gaus");
	tau->Draw();
	gPad->Update();
    
    C6->SaveAs(Form("%s/taus_E.png", outdir));
    C6->SaveAs(Form("%s/taus_E.pdf", outdir));
    
    // Canvas 7: 
    TCanvas* C7 = new TCanvas("C7", "pdg", 1200, 500);
    C7->SetGrid();
    C7->Divide(2,1);
    C7->cd(1);
    
    pdg1->SetLineColor(kAzure);
    pdg1->SetLineWidth(2);
    pdg1->SetTitle("Reco 1 #pi + 1 #gamma NON #rho events : MCparents PDG;PDG;counts");
    pdg1->Draw("HIST");
    
    C7->cd(2);
    
    pdg2->SetLineColor(kOrange+10);
    pdg2->SetLineWidth(2);
    pdg2->SetTitle("Reco 1 #pi + 2 #gamma NON #rho events : MCparents PDG;PDG;counts");
    pdg2->Draw("HIST");
    
    
    C7->SaveAs(Form("%s/pdg.png", outdir));
    C7->SaveAs(Form("%s/pdg.pdf", outdir));
    
    
    
    // Clean up 
    // f->Close();
    // delete f;
}


