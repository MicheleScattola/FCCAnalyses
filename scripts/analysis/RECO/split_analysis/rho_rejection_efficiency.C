#include <iostream>
#include "TFile.h"
#include "TH1I.h"
#include <vector>
#include <iomanip>
#include "ROOT/RDataFrame.hxx"

using namespace ROOT;

void rho_rejection_efficiency() {
    
    // =========================================================================
    // 1. CONFIGURATION
    // =========================================================================
    const std::string input_file = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/split_analysis/output_merged.root";
    const std::string treeName = "events";
    
    std::cout << ">>> Opening file: " << input_file << std::endl;
    
    // =========================================================================
    // 2. CREATE DATAFRAME AND BUILD HISTOGRAMS
    // =========================================================================
    RDataFrame df(treeName, input_file);
    
    auto h_reject_1 = df.Histo1D({"h_rho_mass_reject_1", "rho_mass_reject_1", 6, -0.5, 5.5}, "rho_mass_reject_1");
    auto h_reject_2 = df.Histo1D({"h_rho_mass_reject_2", "rho_mass_reject_2", 6, -0.5, 5.5}, "rho_mass_reject_2");
    
    auto h_angle_1 = df.Histo1D({"h_rho_angle_reject_1", "rho_angle_reject_1", 6, -0.5, 5.5}, "rho_angle_reject_1");
    auto h_angle_2 = df.Histo1D({"h_rho_angle_reject_2", "rho_angle_reject_2", 6, -0.5, 5.5}, "rho_angle_reject_2");
    
    std::cout << ">>> Histograms built via DataFrame" << std::endl;
    
    // =========================================================================
    // 3. MAP MC TYPE VALUES TO CATEGORY NAMES
    // =========================================================================
    std::vector<std::string> categories = {"other", "muons", "electrons", "pions", "rho", "a1"};
    
    // =========================================================================
    // 4. EXTRACT BIN ENTRIES
    // =========================================================================
    std::vector<int> bins_1(6), bins_2(6);
    
    for (int b = 0; b <= 5; ++b) {
        bins_1[b] = static_cast<int>(h_reject_1->GetBinContent(b + 1));  // bin+1 because ROOT bins are 1-indexed
        bins_2[b] = static_cast<int>(h_reject_2->GetBinContent(b + 1));
    }
    
    // =========================================================================
    // 6. PRINT RESULTS FOR rho_mass_reject_1
    // =========================================================================
    std::cout << "\n========================================" << std::endl;
    std::cout << "rho_mass_reject_1 Histogram" << std::endl;
    std::cout << "========================================" << std::endl;
    
    int total_1 = 0;
    int background_1 = 0;
    
    std::cout << std::setw(6) << "Bin"
              << std::setw(15) << "Category"
              << std::setw(15) << "Entries" << std::endl;
    std::cout << std::string(36, '-') << std::endl;
    
    for (int b = 0; b <= 5; ++b) {
        int count = bins_1[b];
        total_1 += count;
        if (b != 4) background_1 += count;  // b=4 is rho
        std::cout << std::setw(6) << b
                  << std::setw(15) << categories[b]
                  << std::setw(15) << count << std::endl;
    }
    std::cout << std::string(36, '-') << std::endl;
    
    double eff_1 = (total_1 > 0) ? (double)background_1 / total_1 : 0.0;
    std::cout << "Total entries: " << total_1 << std::endl;
    std::cout << "Background (non-rho): " << background_1 << std::endl;
    std::cout << "Rejection Efficiency: " << std::fixed << std::setprecision(6) << eff_1 << std::endl;
    
    // =========================================================================
    // 7. PRINT RESULTS FOR rho_mass_reject_2
    // =========================================================================
    std::cout << "\n========================================" << std::endl;
    std::cout << "rho_mass_reject_2 Histogram" << std::endl;
    std::cout << "========================================" << std::endl;
    
    int total_2 = 0;
    int background_2 = 0;
    
    std::cout << std::setw(6) << "Bin"
              << std::setw(15) << "Category"
              << std::setw(15) << "Entries" << std::endl;
    std::cout << std::string(36, '-') << std::endl;
    
    for (int b = 0; b <= 5; ++b) {
        int count = bins_2[b];
        total_2 += count;
        if (b != 4) background_2 += count;  // b=4 is rho
        std::cout << std::setw(6) << b
                  << std::setw(15) << categories[b]
                  << std::setw(15) << count << std::endl;
    }
    std::cout << std::string(36, '-') << std::endl;
    
    double eff_2 = (total_2 > 0) ? (double)background_2 / total_2 : 0.0;
    std::cout << "Total entries: " << total_2 << std::endl;
    std::cout << "Background (non-rho): " << background_2 << std::endl;
    std::cout << "Rejection Efficiency: " << std::fixed << std::setprecision(6) << eff_2 << std::endl;
    
    
    // =========================================================================
    // 8. EXTRACT BIN ENTRIES FOR ANGLE REJECTION
    // =========================================================================
    std::vector<int> angle_bins_1(6), angle_bins_2(6);
    
    for (int b = 0; b <= 5; ++b) {
        angle_bins_1[b] = static_cast<int>(h_angle_1->GetBinContent(b + 1));
        angle_bins_2[b] = static_cast<int>(h_angle_2->GetBinContent(b + 1));
    }
    
    // =========================================================================
    // 9. PRINT RESULTS FOR rho_angle_reject_1
    // =========================================================================
    std::cout << "\n========================================" << std::endl;
    std::cout << "rho_angle_reject_1 Histogram" << std::endl;
    std::cout << "========================================" << std::endl;
    
    int total_angle_1 = 0;
    int background_angle_1 = 0;
    
    std::cout << std::setw(6) << "Bin"
              << std::setw(15) << "Category"
              << std::setw(15) << "Entries" << std::endl;
    std::cout << std::string(36, '-') << std::endl;
    
    for (int b = 0; b <= 5; ++b) {
        int count = angle_bins_1[b];
        total_angle_1 += count;
        if (b != 4) background_angle_1 += count;  // b=4 is rho
        std::cout << std::setw(6) << b
                  << std::setw(15) << categories[b]
                  << std::setw(15) << count << std::endl;
    }
    std::cout << std::string(36, '-') << std::endl;
    
    double eff_angle_1 = (total_angle_1 > 0) ? (double)background_angle_1 / total_angle_1 : 0.0;
    std::cout << "Total entries: " << total_angle_1 << std::endl;
    std::cout << "Background (non-rho): " << background_angle_1 << std::endl;
    std::cout << "Rejection Efficiency: " << std::fixed << std::setprecision(6) << eff_angle_1 << std::endl;
    
    // =========================================================================
    // 10. PRINT RESULTS FOR rho_angle_reject_2
    // =========================================================================
    std::cout << "\n========================================" << std::endl;
    std::cout << "rho_angle_reject_2 Histogram" << std::endl;
    std::cout << "========================================" << std::endl;
    
    int total_angle_2 = 0;
    int background_angle_2 = 0;
    
    std::cout << std::setw(6) << "Bin"
              << std::setw(15) << "Category"
              << std::setw(15) << "Entries" << std::endl;
    std::cout << std::string(36, '-') << std::endl;
    
    for (int b = 0; b <= 5; ++b) {
        int count = angle_bins_2[b];
        total_angle_2 += count;
        if (b != 4) background_angle_2 += count;  // b=4 is rho
        std::cout << std::setw(6) << b
                  << std::setw(15) << categories[b]
                  << std::setw(15) << count << std::endl;
    }
    std::cout << std::string(36, '-') << std::endl;
    
    double eff_angle_2 = (total_angle_2 > 0) ? (double)background_angle_2 / total_angle_2 : 0.0;
    std::cout << "Total entries: " << total_angle_2 << std::endl;
    std::cout << "Background (non-rho): " << background_angle_2 << std::endl;
    std::cout << "Rejection Efficiency: " << std::fixed << std::setprecision(6) << eff_angle_2 << std::endl;


    // =========================================================================
    // 11. SUMMARY
    // =========================================================================
    std::cout << "\n========================================" << std::endl;
    std::cout << "SUMMARY" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "rho_mass_reject_1: Background/Total = " << background_1 << "/" << total_1 
              << " = " << std::fixed << std::setprecision(6) << eff_1 << std::endl;
    std::cout << "rho_mass_reject_2: Background/Total = " << background_2 << "/" << total_2 
              << " = " << std::fixed << std::setprecision(6) << eff_2 << std::endl;
    std::cout << "rho_angle_reject_1: Background/Total = " << background_angle_1 << "/" << total_angle_1 
              << " = " << std::fixed << std::setprecision(6) << eff_angle_1 << std::endl;
    std::cout << "rho_angle_reject_2: Background/Total = " << background_angle_2 << "/" << total_angle_2 
              << " = " << std::fixed << std::setprecision(6) << eff_angle_2 << std::endl;
    
    
    std::cout << "\n>>> Analysis completed." << std::endl;
}
