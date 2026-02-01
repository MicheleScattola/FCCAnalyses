// -*- C++ -*-
//
/** FCCAnalysis module: Z -> tau tau events
 *
 * \file Ztautau.cc
 * \author Michele Scattola <michele.scattola@studenti.unimi.it>
 *
 * Description:
 *   [Analysis header for Z->tau tau events]
 */

#ifndef Ztautau_Ztautau_h
#define Ztautau_Ztautau_h

#include "FCCAnalyses/Algorithms.h"
#include "FCCAnalyses/MCParticle.h"
#include "FCCAnalyses/ReconstructedParticle.h"
#include "FCCAnalyses/VertexingUtils.h"
#include "Math/Vector4D.h"
#include "ROOT/RVec.hxx"
#include "TLorentzVector.h"
#include "edm4hep/ReconstructedParticle.h"

#include <iostream>

using namespace FCCAnalyses;
using namespace ROOT::VecOps;
using namespace ROOT::Math;

namespace Ztautau {

namespace rv = ROOT::VecOps;

const double SM_TAU = 1.77686; // tau mass in GeV
const double SM_PI = 0.13957039;
const double SM_sin2thetaW = 0.2315;
const double gv_ga = 1 - 4 * SM_sin2thetaW;
const double SM_Atau = 2 * gv_ga / (1 + gv_ga * gv_ga);
const double SQRTS = 91.188;    // Z pole energy
const double E_TAU = SQRTS / 2; // tau energy at Z pole

//===================================
// myEVENT
//===================================
struct myEvent {

  // RECO
  int n_mu = 0, n_el = 0, n_pi = 0, n_ph = 0; // particle counts
  double m_RecoCharge = 0.;                   // total charge in hemisphere
  double m_RecoEnergy = 0.;                   // total energy in hemisphere
  double m_RecoMass = 0.; // total invariant mass in hemisphere
  RVec<TLorentzVector> m_muP4, m_elP4, m_piP4,
      m_phP4;            // particle TLorentzVectors
  int m_type = 0;        // event reco type
  int m_type_before = 0; // event reco type before optimal variable constraints
  double m_omega = 0;    // optimal variable binned

  // MC
  int mc_tau_index = -1;        // tau MC index
  int mc_type = 0;              // event mc type
  double mc_Ptau = 0.;          // tau polarization from MC
  double mc_weight_plus = 1.0;  // reweighting for h = +1
  double mc_weight_minus = 1.0; // reweighting for h = -1
  double mc_omega = 0.;         // optimal variable for each channel
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

  // debug
  int m_debug = 0;
  int m_debug_mass = 0;
  double thrust_costheta = -999;
  double thrust_phi = -999;
  double thrust_costheta_hemi = -999;

  bool is_dressed = false;
};

//===================================
// EVENT STRUCT
//===================================

// return event struct
RVec<myEvent> myget_event(const RVec<int> &mu_ids, const RVec<int> &el_ids,
                          const RVec<int> &pi_ids, const RVec<int> &ph_ids,
                          const RVec<edm4hep::ReconstructedParticleData> &rps,
                          const RVec<double> &rps_costheta,
                          const RVec<edm4hep::MCParticleData> &mc,
                          const RVec<int> &daughters,
                          const RVec<int> &rp2mc_idx,
                          const RVec<double> &thrust_vector);

// return signed thrust costheta per hemisphere
RVec<double> get_thrustcostheta_hemi(const RVec<myEvent> &evs);

// ==========================================
// helper to fill P4 and add energy & charge
void fill_collection(myEvent &ev, const auto &input_particles,
                     RVec<TLorentzVector> &out_p4);

// ==========================================
// helper to retrieve indices based on hemi sign
RVec<int> get_idx(const bool &hemi, const RVec<double> &costheta,
                  const RVec<int> &ids);

//===================================
// custom getThrustPointing using charge instead of energy
RVec<double> getThrustPointing(const RVec<double> &charge,
                               const RVec<double> &thrust,
                               const RVec<double> &costheta);

//===================================
// select pion indices as hadron with mass = 0.13957039
RVec<int> sel_pions_id(const RVec<edm4hep::ReconstructedParticleData> &in,
                       const int charge);

//===================================
// select elements given indices
template <typename T>
inline RVec<T> get_elements_by_index(const RVec<T> &A, const RVec<int> &B) {
  RVec<T> out;
  out.reserve(B.size());
  const int N = static_cast<int>(A.size());
  for (int idx : B) {
    if (idx >= 0 && idx < N)
      out.push_back(A[idx]);
  }
  return out;
}

//===================================
// EVENT CLASSIFICATION
//===================================

// ==========================================
int classify_lep(myEvent &ev, const RVec<int> &mu_idx, const RVec<int> &el_idx,
                 const RVec<int> &rp2mc_idx);

int classify_pion(myEvent &ev, const RVec<int> &pi_idx,
                  const RVec<int> &rp2mc_idx);

int classify_MC(const RVec<int> &pdgs);

// ==========================================
RVec<int> get_type_safe(const RVec<myEvent> &evs);

// ==========================================
// RECO OPTIMAL VARIABLE CALCULATION
// ==========================================

void reco_omega(myEvent &ev);

double calculate_omega_rho(myEvent &ev, const TLorentzVector &p4_tau,
                           const TLorentzVector &p4_rho,
                           const TLorentzVector &p4_pip,
                           const TLorentzVector &p4_pi0);

// ==========================================
// EXTRACT VARIABLES
// ==========================================

RVec<double> get_optimal(RVec<myEvent> &evs, const int mc_type,
                         const bool bool_mc, const int reco_type,
                         const bool bool_reco);

// =========================================
RVec<double> get_reco_x(RVec<myEvent> &evs, const int mc_type,
                        const bool bool_mc, const int reco_type,
                        const bool bool_reco);

// ==========================================
RVec<double> get_weights(const int sign, const RVec<myEvent> &evs,
                         const int mc_type, const bool bool_mc,
                         const int reco_type, const bool bool_reco);

// ==========================================
// RE-WEIGHTING FUNCTIONS
// ==========================================
double calc_Ptau(const TLorentzVector &p4_tau);

double GetCosThetaStar(const TLorentzVector &p4_tau_lab,
                       const TLorentzVector &p4_pi_lab);

void lepton_weight(myEvent &ev, const RVec<edm4hep::MCParticleData> &mc,
                   const RVec<int> &daughters);

void pion_weight(myEvent &ev, const RVec<edm4hep::MCParticleData> &mc,
                 const RVec<int> &daughters);

void rho_weight(myEvent &ev, const RVec<edm4hep::MCParticleData> &mc,
                const RVec<int> &daughters);

void a1_weight(myEvent &ev, const RVec<edm4hep::MCParticleData> &mc,
               const RVec<int> &daughters);

void new_rho_weight(myEvent &ev, const RVec<edm4hep::MCParticleData> &mc,
                    const RVec<int> &daughters);

// omega_rho with p4 angles instead of kinematic variables, possible only for MC
double geometric_omega_rho(myEvent &ev, TLorentzVector &p4_tau,
                           const TLorentzVector &p4_rho,
                           const TLorentzVector &p4_pip,
                           const TLorentzVector &p4_pi0);

// ==========================================
// additional
// ==========================================

// ==========================================
RVec<double> get_invariant_mass(const RVec<myEvent> &evs, const int mc_type,
                                const bool bool_mc, const int reco_type,
                                const bool bool_reco);

// ==========================================
RVec<double> get_mass_pull(const RVec<myEvent> &evs, const int mc_type,
                           const bool bool_mc, const int reco_type,
                           const bool bool_reco);

// ==========================================
RVec<double> get_Ptau(const RVec<myEvent> &evs, const int mc_type,
                      const bool bool_mc, const int reco_type,
                      const bool bool_reco, const bool masscheck,
                      const bool asymmetric);

// ==========================================
RVec<double> get_omega_rho(RVec<myEvent> &evs, const int mc_type,
                           const bool bool_mc, const int reco_type,
                           const bool bool_reco, const bool masscheck,
                           const bool asymmetric);

// ==========================================
// DEBUG
// ==========================================

RVec<int> get_debug_mass(const RVec<myEvent> &evs, const int reco_type,
                         const int n_photons, const bool bool_ph);

RVec<int> get_debug(const RVec<myEvent> &evs, const int reco_type,
                    const int n_photons);

// =========================================
// NOT IN USE
// =========================================
RVec<double> get_lepton_x(const RVec<myEvent> &evs, const int mc_type,
                          const bool bool_mc, const int reco_type,
                          const bool bool_reco, const bool masscheck,
                          const bool asymmetric);

RVec<double> get_hadron_e(const RVec<myEvent> &evs, const int mc_type,
                          const bool bool_mc, const int reco_type,
                          const bool masscheck);

// ==========================================
RVec<double> get_photon_e(const RVec<myEvent> &evs, const int mc_type,
                          const bool bool_mc, const int reco_type,
                          const bool bool_reco, const bool masscheck,
                          const bool asymmetric);
// ==========================================
RVec<double> get_dressed_e(const RVec<myEvent> &evs, const int mc_type,
                           const bool bool_mc, const int reco_type,
                           const bool bool_reco, const bool masscheck,
                           const bool asymmetric);
// ==========================================
RVec<double> get_rp2mc_e(const RVec<myEvent> &evs, const int mc_type,
                         const bool bool_mc, const int reco_type,
                         const bool bool_reco, const bool masscheck,
                         const bool asymmetric);

// ==========================================
RVec<double> get_MCdaughter_e(const RVec<myEvent> &evs, const int mc_type,
                              const bool bool_mc, const int reco_type,
                              const bool bool_reco, const bool masscheck,
                              const bool asymmetric);
// ==========================================
RVec<double> get_MCdaughter_x(const RVec<myEvent> &evs, const int mc_type,
                              const bool bool_mc, const int reco_type,
                              const bool bool_reco, const bool masscheck,
                              const bool asymmetric);
// ==========================================
RVec<double> get_MCdaughter_mass(const RVec<myEvent> &evs, const int mc_type,
                                 const bool bool_mc, const int reco_type,
                                 const bool bool_reco);

// ==========================================
// get reco optimal in specifiec costheta min and max
RVec<double> get_reco_x_theta(RVec<myEvent> &evs, const int mc_type,
                              const bool bool_mc, const int reco_type,
                              const bool bool_reco, const double costheta_min,
                              const double costheta_max);

// ==========================================
// get weights in specific costheta min and max
RVec<double> get_weights_theta(const int sign, const RVec<myEvent> &evs,
                               const int mc_type, const bool bool_mc,
                               const int reco_type, const bool bool_reco,
                               const double costheta_min,
                               const double costheta_max);

} // namespace Ztautau

#endif
