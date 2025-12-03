{
    gROOT->Reset();

    // Input and output paths (edit if needed)
    const char* infile = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/histmaker/p8_ee_Ztautau_ecm91.root";
    const char* outdir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/plots/single_pi";

    TFile* f = TFile::Open(infile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: could not open file: " << infile << std::endl;
        return;
    }

    // Retrieve histograms 
    TH1F* ea = (TH1F*)f->Get("el1g_all");
    TH1F* ma = (TH1F*)f->Get("mu1g_all");
    TH1F* ta = (TH1F*)f->Get("tau1g_all");
    TH1F* pa = (TH1F*)f->Get("ph1g_all");
    TH1F* pi0a = (TH1F*)f->Get("pi01g_all");
    
    TH1F* et = (TH1F*)f->Get("el1g_true");
    TH1F* mt = (TH1F*)f->Get("mu1g_true");
    TH1F* tt = (TH1F*)f->Get("tau1g_true");
    TH1F* pt = (TH1F*)f->Get("ph1g_true");
    TH1F* pi0t = (TH1F*)f->Get("pi01g_true");
    
    // Check histograms existence (all 1γ ALL)
    if (!ea)   std::cerr << "ERROR: histogram 'el1g_all' not found\n";
    if (!ma)   std::cerr << "ERROR: histogram 'mu1g_all' not found\n";
    if (!ta)   std::cerr << "ERROR: histogram 'tau1g_all' not found\n";
    if (!pa)   std::cerr << "ERROR: histogram 'ph1g_all' not found\n";
    if (!pi0a) std::cerr << "ERROR: histogram 'pi01g_all' not found\n";

    // Check histograms existence (1γ TRUE)
    if (!et)   std::cerr << "ERROR: histogram 'el1g_true' not found\n";
    if (!mt)   std::cerr << "ERROR: histogram 'mu1g_true' not found\n";
    if (!tt)   std::cerr << "ERROR: histogram 'tau1g_true' not found\n";
    if (!pt)   std::cerr << "ERROR: histogram 'ph1g_true' not found\n";
    if (!pi0t) std::cerr << "ERROR: histogram 'pi01g_true' not found\n";
    


    // Create output directory if missing (system call)
    gSystem->Exec(Form("mkdir -p %s", outdir));

    // -------------------------
    // Single canvas with all plots
    TCanvas* C = new TCanvas("C", "all 1 gamma", 800, 500);
    C->SetGrid();
    C->cd();

    pi0a->SetTitle("Reco 1#pi + 1#gamma ALL events: MCparentes PDG contributions;Energy_{#gamma} [GeV];Entries");
    pi0a->SetLineColor(kRed);
    pi0a->SetLineWidth(2);
    pi0a->Scale(1/10.);
    pi0a->Draw("HIST");
    
    ea->SetLineColor(kGreen);
    ea->SetLineWidth(2);
    ea->Draw("HIST SAME");
    
    ma->SetLineColor(kOrange);
    ma->SetLineWidth(2);
    ma->Draw("HIST SAME");
    
    ta->SetLineColor(kBlue);
    ta->SetLineWidth(2);
    ta->Draw("HIST SAME");
    
    pa->SetLineColor(kCyan);
    pa->SetLineWidth(2);
    pa->Draw("HIST SAME");
    
    // Legend
    TLegend* leg = new TLegend(0.70,0.55,0.88,0.88);
    leg->AddEntry(ea, "e", "l");
    leg->AddEntry(ma, "#mu", "l");
    leg->AddEntry(ta, "#tau", "l");
    leg->AddEntry(pt, "#gamma", "l");
    leg->AddEntry(pi0a, "#pi^{0} (scaled /10)", "l");
    leg->Draw();

    C->Update();
    C->SaveAs(Form("%s/all_1gamma.png", outdir));
    C->SaveAs(Form("%s/all_1gamma.pdf", outdir));
    
    // -------------------------
    // Single canvas with all plots
    TCanvas* C1 = new TCanvas("C1", "true tau to pi", 800, 500);
    C1->SetGrid();
    C1->cd();

    tt->SetTitle("Reco 1#pi + 1#gamma TrueMC events: MCparents PDG contributions;Energy_{#gamma} [GeV];Entries");
    tt->SetLineColor(kBlue);
    tt->SetLineWidth(2);
    tt->Scale(1/2.);
    tt->Draw("HIST");
    
    et->SetLineColor(kGreen);
    et->SetLineWidth(2);
    et->Draw("HIST SAME");
    
    mt->SetLineColor(kOrange);
    mt->SetLineWidth(2);
    mt->Draw("HIST SAME");
    
    pt->SetLineColor(kCyan);
    pt->SetLineWidth(2);
    pt->Draw("HIST SAME");
    
    pi0t->SetLineColor(kRed);
    pi0t->SetLineWidth(2);
    pi0t->Draw("HIST SAME");
    
    // Legend
    TLegend* leg1 = new TLegend(0.70,0.55,0.88,0.88);
    leg1->AddEntry(et, "e", "l");
    leg1->AddEntry(mt, "#mu", "l");
    leg1->AddEntry(tt, "#tau (scaled /2)", "l");
    leg1->AddEntry(pt, "#gamma", "l");
    leg1->AddEntry(pi0t, "#pi^{0} ", "l");
    leg1->Draw();

    C1->Update();
    C1->SaveAs(Form("%s/true_1gamma.png", outdir));
    C1->SaveAs(Form("%s/true_1gamma.pdf", outdir));
    
    
    
    // Clean up 
    // f->Close();
    // delete f;
}


