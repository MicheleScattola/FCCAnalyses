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

const float SM_TAU = 1.77686; // tau mass in GeV
const float SM_PI_CHARGED = 0.13957039;
const float SM_sin2thetaW = 0.2312;
const float gv_ga = -1 + 4 * SM_sin2thetaW;
const float SM_Atau = 2 * gv_ga / (1 + gv_ga * gv_ga);
const float SM_P_TAU =
    -0.1421;              // tau polarization in Z decays at sqrt(s) = 91.2 GeV
const float SQRTS = 91.2; // Z pole energy
const float E_TAU = SQRTS / 2; // tau energy at Z pole
//===================================
// custom getThrustPointing using charge instead of energy
RVec<float> getThrustPointing(const RVec<float> &charge,
                              const RVec<float> &thrust,
                              const RVec<float> &costheta);

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
  float m_RecoCharge = 0.;                    // total charge in hemisphere
  float m_RecoEnergy = 0.;                    // total energy in hemisphere
  float m_RecoMass = 0.;                      // total invariant mass in hemisphere
  RVec<TLorentzVector> m_muP4, m_elP4, m_piP4,
      m_phP4;     // particle TLorentzVectors
  int m_type = 0; // event reco type

  // MC
  int m_tauMCindex = -1;        // tau MC index
  int m_MCtype = 0;             // event mc type
  float m_MCweight_plus = 1.0;  // reweighting for h = +1
  float m_MCweight_minus = 1.0; // reweighting for h = -1
  bool m_found = false;
  RVec<TLorentzVector> mc_tauP4;
  RVec<TLorentzVector> mc_piP4;

  // debug
  int m_debug = 0; 
  float m_invariant_mass = 0.0;
};

// return event struct
RVec<myEvent> myget_event(const RVec<int> &mu_ids, const RVec<int> &el_ids,
                          const RVec<int> &pi_ids, const RVec<int> &ph_ids,
                          const RVec<edm4hep::ReconstructedParticleData> &rps,
                          const RVec<float> &rps_costheta,
                          const RVec<edm4hep::MCParticleData> &mc,
                          const RVec<int> &daughters);
// ==========================================
// helper to fill P4 and add energy & charge
auto fill_collection = (const auto &input_particles,
                            RVec<TLorentzVector> &out_p4);
// ==========================================
int classify_lep(const myEvent &ev);
int classify_pion(
    myEvent &ev); // NOT CONST, need to add invariant mass information
int classify_MC(const RVec<int> &pdgs);
// ==========================================
RVec<int> get_type_safe(const RVec<myEvent> &evs);
// ==========================================
// RE-WEIGHTING FUNCTIONS
// ==========================================
float calc_Ptau(const TLorentzVector &p4_tau);

float GetCosThetaStar(const TLorentzVector &p4_tau_lab,
                      const TLorentzVector &p4_pi_lab);

void pion_weight(myEvent &ev, const RVec<edm4hep::MCParticleData> &mc,
                 const RVec<int> &daughters);

void rho_weight(myEvent &ev, const RVec<edm4hep::MCParticleData> &mc,
                const RVec<int> &daughters);

void a1_weight(myEvent &ev, const RVec<edm4hep::MCParticleData> &mc,
               const RVec<int> &daughters);

// ==========================================
// EXTRACT OPTIMAL VARIABLES
// ==========================================
RVec<float> get_lepton_e(const RVec<myEvent> &evs, const int mc_type,
                         const bool bool_mc, const int reco_type,
                         const bool bool_reco);

RVec<float> get_hadron_e(const RVec<myEvent> &evs, const int mc_type,
                         const bool bool_mc, const int reco_type,
                         const bool bool_reco);
// ==========================================
RVec<float> get_photon_e(const RVec<myEvent> &evs, const int mc_type,
                         const bool bool_mc, const int reco_type,
                         const bool bool_reco);
// ==========================================
RVec<float> get_MCpi_e(const RVec<myEvent> &evs, const int mc_type,
                       const bool bool_mc, const int reco_type,
                       const bool bool_reco);
// ==========================================
RVec<float> get_weights(const int sign, const RVec<myEvent> &evs,
                        const int mc_type, const bool bool_mc,
                        const int reco_type, const bool bool_reco);

// ==========================================
// MASKS AND FILTERS
// ==========================================
RVec<float> get_invariant_mass(const RVec<myEvent> &evs, const int mc_type,
                               const bool bool_mc, const int reco_type,
                               const bool bool_reco);
RVec<int> get_pi_mask(const RVec<myEvent> &evs);
RVec<int> get_weight_mask(const RVec<myEvent> &evs);

} // namespace Ztautau

#endif
