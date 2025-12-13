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
const double SM_PI_CHARGED = 0.13957039;
const double SM_sin2thetaW = 0.23148;
const double gv_ga = 1 - 4 * SM_sin2thetaW;
const double SM_Atau = 2 * gv_ga / (1 + gv_ga * gv_ga);
const double SQRTS = 91.2; // Z pole energy
const double E_TAU = SQRTS / 2; // tau energy at Z pole
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
// myEVENT
//===================================
struct myEvent {

  // RECO
  int n_mu = 0, n_el = 0, n_pi = 0, n_ph = 0; // particle counts
  double m_RecoCharge = 0.;                    // total charge in hemisphere
  double m_RecoEnergy = 0.;                    // total energy in hemisphere
  double m_RecoMass = 0.;                      // total invariant mass in hemisphere
  RVec<TLorentzVector> m_muP4, m_elP4, m_piP4,
      m_phP4;     // particle TLorentzVectors
  int m_type = 0; // event reco type

  // MC
  int m_tauMCindex = -1;        // tau MC index
  int m_MCtype = 0;             // event mc type
  double m_MCPtau = 0.;        // tau polarization from MC
  double m_MCweight_plus = 1.0;  // reweighting for h = +1
  double m_MCweight_minus = 1.0; // reweighting for h = -1
  bool m_found = false;
  TLorentzVector mc_tauP4;    // P4 of mc tau
  TLorentzVector mc_daughterP4;   // P4 of the mc daughter (whole resonance in case)
  RVec<int> mc_daughters;   // vector of daughters pdgs
  int mc_RP2MC_id = -1;   // corresponding mc id for the 1prong reco particle
  double mc_RP2MC_e = 0.;   // correponding energy
  double mc_daughterMass = 0.;    // mc invariant mass of the mc daughter

  // debug
  int m_debug = 0; 
  int m_debug_mass = 0;
  double thrust_costheta = -999;
  double thrust_phi = -999;
};

// return event struct
RVec<myEvent> myget_event(const RVec<int> &mu_ids, const RVec<int> &el_ids,
                          const RVec<int> &pi_ids, const RVec<int> &ph_ids,
                          const RVec<edm4hep::ReconstructedParticleData> &rps,
                          const RVec<double> &rps_costheta,
                          const RVec<edm4hep::MCParticleData> &mc,
                          const RVec<int> &daughters,
                          const RVec<int> &rp2mc_idx,
                          const double &thrust_costheta,
                          const double &thrust_phi);
// ==========================================
// helper to fill P4 and add energy & charge
void fill_collection (myEvent &ev, const auto &input_particles,
                            RVec<TLorentzVector> &out_p4);
// ==========================================
// helper
RVec<int> get_idx(const bool &hemi, const RVec<double> &costheta, const RVec<int> &ids);
// ==========================================
int classify_lep(myEvent &ev, const RVec<int> &mu_idx, const RVec<int> &el_idx, const RVec<int> &rp2mc_idx);
int classify_pion(myEvent &ev, const RVec<int> &pi_idx, const RVec<int> &rp2mc_idx);
int classify_MC(const RVec<int> &pdgs);
// ==========================================
RVec<int> get_type_safe(const RVec<myEvent> &evs);
// ==========================================
// RE-WEIGHTING FUNCTIONS
// ==========================================
double calc_Ptau(const TLorentzVector &p4_tau);

double GetCosThetaStar(const TLorentzVector &p4_tau_lab,
                      const TLorentzVector &p4_pi_lab);

void pion_weight(myEvent &ev, const RVec<edm4hep::MCParticleData> &mc,
                 const RVec<int> &daughters);

void rho_weight(myEvent &ev, const RVec<edm4hep::MCParticleData> &mc,
                const RVec<int> &daughters);

void a1_weight(myEvent &ev, const RVec<edm4hep::MCParticleData> &mc,
               const RVec<int> &daughters);

// ==========================================
// EXTRACT VARIABLES
// ==========================================
RVec<double> get_lepton_e(const RVec<myEvent> &evs, const int mc_type,
                         const bool bool_mc, const int reco_type,
                         const bool bool_reco, const bool masscheck);

RVec<double> get_hadron_e(const RVec<myEvent> &evs, const int mc_type,
                         const bool bool_mc, const int reco_type,
                         const bool bool_reco, const bool masscheck);
// ==========================================
RVec<double> get_photon_e(const RVec<myEvent> &evs, const int mc_type,
                         const bool bool_mc, const int reco_type,
                         const bool bool_reco, const bool masscheck);
// ==========================================
RVec<double> get_dressed_e(const RVec<myEvent> &evs, const int mc_type,
                         const bool bool_mc, const int reco_type,
                         const bool bool_reco, const bool masscheck);
// ==========================================
RVec<double> get_rp2mc_e(const RVec<myEvent> &evs, const int mc_type,
                         const bool bool_mc, const int reco_type,
                         const bool bool_reco, const bool masscheck);
// ==========================================
RVec<double> get_weights(const int sign, const RVec<myEvent> &evs,
                        const int mc_type, const bool bool_mc,
                        const int reco_type, const bool bool_reco, const bool masscheck); 
// ==========================================          
RVec<double> get_invariant_mass(const RVec<myEvent> &evs, const int mc_type,
                               const bool bool_mc, const int reco_type,
                               const bool bool_reco);
// ==========================================
RVec<double> get_mass_pull(const RVec<myEvent> &evs, const int mc_type,
                       const bool bool_mc, const int reco_type,
                       const bool bool_reco);
// ==========================================
RVec<double> get_MCdaughter_e(const RVec<myEvent> &evs, const int mc_type,
                       const bool bool_mc, const int reco_type,
                       const bool bool_reco, const bool masscheck);
// ==========================================
RVec<double> get_MCdaughter_x(const RVec<myEvent> &evs, const int mc_type,
                       const bool bool_mc, const int reco_type,
                       const bool bool_reco, const bool masscheck);
// ==========================================
RVec<double> get_MCdaughter_mass(const RVec<myEvent> &evs, const int mc_type,
                       const bool bool_mc, const int reco_type,
                       const bool bool_reco);
// ==========================================
// MASKS AND FILTERS
// ==========================================
RVec<int> get_pi_mask(const RVec<myEvent> &evs);
RVec<int> get_weight_mask(const RVec<myEvent> &evs);
RVec<int> get_debug(const RVec<myEvent> &evs,
                               const int mc_type, const bool bool_mc,
                               const int reco_type, const bool bool_reco,
                               const int debug_mass, const bool bool_debug_mass);

RVec<RVec<int>> get_debug_daughters(const RVec<myEvent> &evs,
                               const int mc_type, const bool bool_mc,
                               const int reco_type, const bool bool_reco);
RVec<double> get_debug_costheta(const RVec<myEvent> &evs,
                               const int mc_type, const bool bool_mc,
                               const int reco_type, const bool bool_reco);
RVec<double> get_debug_phi(const RVec<myEvent> &evs,
                               const int mc_type, const bool bool_mc,
                               const int reco_type, const bool bool_reco);
} // namespace Ztautau

#endif
