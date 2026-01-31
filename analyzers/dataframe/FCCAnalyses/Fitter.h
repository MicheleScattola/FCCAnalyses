#ifndef FITTER_H
#define FITTER_H

#include <iostream>
#include <string>
#include <vector>

namespace Fitter {

// Result struct
struct myFit {
  double P_tau;   // Polarization
  double P_err;   // Error
  double f_plus;  // Fraction H+ (Template only)
  double f_minus; // Fraction H- (Template only)
  bool success;
  std::string method; // "Template" or "Analytic"
};

// =========================================================
// OLD METHOD: FRACTION FIT (Fit N+, N- then Error Prop)
// =========================================================
myFit fit_fraction(const std::string &infile_data,
                   const std::string &infile_templates,
                   const std::string &outdir,
                   const std::string &output_filename,
                   const std::string &treeName, const std::string &dataColName,
                   const std::string &name_plus, const std::string &name_minus,
                   const std::string &plot_title,
                   const std::string &x_axis_title);

// =========================================================
// TEMPLATE FIT (generic)
// =========================================================
myFit fit(const std::string &infile_data, const std::string &infile_templates,
          const std::string &outdir, const std::string &output_filename,
          const std::string &treeName, const std::string &dataColName,
          const std::string &name_plus, const std::string &name_minus,
          const std::string &plot_title, const std::string &x_axis_title);

// =========================================================
// TEMPLATE FIT WITHOUT PLOTTING (for batch processing)
// =========================================================
myFit fit_no_plot(const std::string &infile_data,
                  const std::string &infile_templates,
                  const std::string &treeName, const std::string &dataColName,
                  const std::string &name_plus, const std::string &name_minus);

// =========================================================
// TEMPLATE FIT WITH DATAFRAME FILTERING (for angular binning)
// =========================================================
myFit fit_filtered(const std::string &infile_data,
                   const std::string &infile_templates,
                   const std::string &treeName, const std::string &dataColName,
                   const std::string &name_plus, const std::string &name_minus);

// =========================================================
// ANALYTIC FIT (leptons)
// =========================================================
myFit fit(const std::string &infile_data, const std::string &outdir,
          const std::string &output_filename, const std::string &treeName,
          const std::string &dataColName, const std::string &plot_title,
          const std::string &x_axis_title, double xmin, double xmax);

} // namespace Fitter

#endif
