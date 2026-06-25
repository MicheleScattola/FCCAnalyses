# Analysis and code implemented
Code was written by Michele Scattola, if you have any questions don't hesitate to write at michele.scattola@studenti.unimi.it
The main headers are found in `FCCAnalyses/analyzer/dataframe/FCCAnalyses/Ztautau.h` and `FCCAnalyses/analyzer/dataframe/src/Ztautau.cc`.

The analysis creates instances of a `struct myEvent`. These hold information about the reconstructed and Montecarlo level event, such as classification, P4 vectors and debug variables.

### struct myEvent
- Reconstructed (RECO) information is stored in the following variables:
  ```c++
  // RECO
  int n_mu = 0, n_el = 0, n_pi = 0, n_ph = 0; // particle counts
  double m_RecoCharge = 0.;                   // total charge in hemisphere
  double m_RecoEnergy = 0.;                   // total energy in hemisphere
  double m_RecoMass = 0.; // total invariant mass in hemisphere
  RVec<TLorentzVector> m_muP4, m_elP4, m_piP4,
      m_phP4;             // particle TLorentzVectors
  TLorentzVector m_tauP4; // reconstructed tau P4
  int m_type = 0;         // event reco type
  int m_type_before = 0;  // event reco type before optimal variable constraints
  double m_omega = 0;     // optimal variable binned
  ```

- Montecarlo (MC) information stores in analogy:
  ```c++
  // MC
  int mc_tau_index = -1;        // tau MC index
  int mc_type = 0;              // event mc type
  double mc_Ptau = 0.;          // tau polarization from MC
  double mc_weight_plus = 1.0;  // reweighting for h = +1
  double mc_weight_minus = 1.0; // reweighting for h = -1
  double mc_omega = 0.;         // optimal variable for each channel
  double mc_charge = 0.;
  bool m_found = false;
  TLorentzVector mc_tauP4; // P4 of mc tau
  TLorentzVector
      mc_daughterP4; // P4 of the mc daughter (whole resonance in case)
  TLorentzVector mc_pi0P4;
  TLorentzVector mc_piP4;
  RVec<int> mc_daughters;     // vector of daughters pdgs
  RVec<int> mc_debug_parents; // vector of debug fathers pdgs
  int mc_RP2MC_id = -1;   // corresponding mc id for the 1prong reco particle
  double mc_RP2MC_e = 0.; // correponding energy
  double mc_daughterMass = 0.; // mc invariant mass of the mc daughter
  int mc_debug_mass = 0;
  ```

- Some debug variables which I found useful during my work:
  ```c++
  // debug
  int m_debug = 0;
  int m_debug_mass = 0;
  double thrust_costheta = -999;
  double thrust_phi = -999;
  double thrust_costheta_hemi = -999;
  TVector3 thrust_vector;

  bool is_dressed = false;
  ```

### Reconstruction
The event is reconstructed in the function `RVec<myEvent> myget_event`. It can be summarized in a few steps.
- build thrust vector (prior step in RDataFrame).
- divide RECO particles in forward and backwards hemispheres w.r.t. thrust.
- implement event identification for each hemisphere for both RECO and MC. (also apply kinematic/invariant mass cuts).
- given the identification calculate the optimal kinematic variable $\omega$ and the histogram re-weighting weight (if needed).
  
### Miscellaneous
- The function `void reco_omega` takes care of calculating the optimal variable, depending on the assigned RECO event identification.
- The function `void rho_weight` showed some unclear behaviour. In order to recover full sensitivity to the polarization the re-weighting variable should be the optimal variable $\omega$, however it produced very instable results. Only partial sensitivity is obtained in the current implementation, where the spin-analyzer is the angle between the rho and the tau, see:
  ```c++
  z = GetCosThetaStar(p4_tau_lab, p4_rho_lab);
  alpha = (SM_TAU * SM_TAU - 2 * mRho * mRho) / (SM_TAU * SM_TAU + 2 * mRho * mRho);
  ev.mc_weight_plus = (1 + alpha * z) / (1 + alpha * Ptau * z);
  ev.mc_weight_minus = (1 - alpha * z) / (1 + alpha * Ptau * z);
  ```
  an optimal approach should implement instead

  ```c++
  ev.mc_weight_plus = = (1 + omega) / (1 + Ptau * omega);
  ev.mc_weight_minus = (1 - omega) / (1 + Ptau * omega);
  ```
  This issue is reflected in the loss of sensitivity to the final measurement, see slide 10 of *presentazione_Scattola.pdf* and chapter 5.3 of *Thesis_Scattola.pdf*