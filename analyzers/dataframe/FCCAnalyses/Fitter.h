#ifndef FITTER_H
#define FITTER_H

#include <string>
#include <vector>
#include <iostream>

namespace Fitter {

    // Result struct
    struct FitResult {
        double P_tau;       // Polarization
        double P_err;       // Error
        double f_plus;      // Fraction H+ (Template only)
        double f_minus;     // Fraction H- (Template only)
        bool success;
        std::string method; // "Template" or "Analytic"
    };

    // =========================================================
    // TEMPLATE FIT (generic)
    // =========================================================
    FitResult fit(const std::string& infile_data,
                  const std::string& infile_templates,
                  const std::string& outdir,
                  const std::string& output_filename, 
                  const std::string& treeName,
                  const std::string& dataColName,     
                  const std::string& name_plus,       
                  const std::string& name_minus,      
                  const std::string& plot_title,      
                  const std::string& x_axis_title     
                  );

    // =========================================================
    // ANALYTIC FIT (leptons)
    // =========================================================
    FitResult fit(const std::string& infile_data,
                  const std::string& outdir,
                  const std::string& output_filename, 
                  const std::string& treeName,
                  const std::string& dataColName,     
                  const std::string& plot_title,      
                  const std::string& x_axis_title,
                  double xmin,
                  double xmax
                  );

} // namespace Fitter

#endif