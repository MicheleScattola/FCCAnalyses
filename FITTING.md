# Fitting methods
Code was written by Michele Scattola, if you have any questions don't hesitate to write at michele.scattola@studenti.unimi.it
The main headers are found in `FCCAnalyses/analyzer/dataframe/FCCAnalyses/Fitter.h` and `FCCAnalyses/analyzer/dataframe/src/Fitter.cc`.
These functions are then called from root macros in order to perform the required fits.

The fitting methods follow directly from the necessity of implementing a histogram re-weighting technique.
I generated templates from Montecarlo data, therefore data is fitted on a linear combination of the polarized templates which then yields the final polarization estimate.

The `struct TemplateFunctor` takes care of accessing the templates and evaluate the current fit iteration from the linear combination of the two. Note that the polarization is obtained via $P = \frac{N_+ - N_-}{N_+ + N_-}$. Given that the templates are normalized (PDFs), there is a normalization to carry which is fixed at the number of events in the histogram. The fit is infact perfored for one parameter only $P$ and not the two values $N_+$ and $N_-$:
``` // Data = Norm * [ ((1+P)/2)*H_plus + ((1-P)/2)*H_minus ] ```

Given the different types of analysis there a couple of fitting functions but they essentially behave in the same way.
- regular generic fit from templates with plotting.
- `fit_no_plot` fitting without plotting, useful for batch processing of many datasets.
- `fit_filtered` is used in angular analysis, where the dataframes are filtered for angular bins (for example in leptonic universality analysis).