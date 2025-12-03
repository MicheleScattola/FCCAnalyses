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
#include "ROOT/RVec.hxx"
#include "TLorentzVector.h"
#include "edm4hep/ReconstructedParticle.h"
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <TStyle.h>

#include <iostream>

using namespace std;
using namespace FCCAnalyses;
using namespace ROOT::VecOps;

namespace Ztautau {

namespace rv = ROOT::VecOps;

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
//===================================
// EVENT CLASSIFICATION
//===================================
//===================================

struct event {
  int m_muN, m_elN, m_piN, m_phN;
  float m_ph_e_max;
  int type;
  RVec<float> m_ph_e;
  float m_sum_ph_e;
  bool m_correct_id = false;
};

// auxiliary function for classification
event classify(int n_mu, int n_el, int n_pi, int n_ph, const RVec<float> &ph_e,
               float ph_e_cutoff);

//===================================
// return event struct
RVec<event> get_event(const RVec<edm4hep::ReconstructedParticleData> &mu,
                      const RVec<float> &mu_costheta,
                      const RVec<edm4hep::ReconstructedParticleData> &el,
                      const RVec<float> &el_costheta,
                      const RVec<edm4hep::ReconstructedParticleData> &pi,
                      const RVec<float> &pi_costheta,
                      const RVec<edm4hep::ReconstructedParticleData> &ph,
                      const RVec<float> &ph_costheta, const float ph_e_cutoff);

// auxiliary function, return RVec of particles energies in passed hemisphere
// sign useful for photons
RVec<int> aux_n(const RVec<float> &legs_costheta);

// auxiliary function, return RVec of particles energies in passed hemisphere
// sign
RVec<float> aux_e(int hemisphere,
                  const RVec<edm4hep::ReconstructedParticleData> &legs,
                  const RVec<float> &legs_costheta);

//===================================
// montecarlo classification
int MC_classified(const RVec<int> &TauMu, const RVec<int> &TauEl,
                  const RVec<int> &TauPi, const RVec<int> &TauRho);

// vector of photon energies given conditions
// bool for TrueMC rho event
// RVec<event> &ev are RecoPart properties of event (per hemisphere)
// int N number of RecoPhotons collected
RVec<float> study_ph(const RVec<int> &MC_event, const RVec<event> &ev, bool foo,
                     int N_ph, bool wrong_events, int type);

// return sum of energies
RVec<float> study_ph_sum(const RVec<int> &MC_event, const RVec<event> &ev,
                         bool foo, int N_ph, bool wrong_events, int type);

//===================================
//===================================
// KINEMATICS
//===================================
//===================================

// building all pairs for pi0 resonance from photons
RVec<edm4hep::ReconstructedParticleData>
pi0_resonance_pairs(const RVec<edm4hep::ReconstructedParticleData> &legs,
                    const float target_mass);

struct pi0_candidate {

  static constexpr float mpi0 = 0.1349768;
  float mass;
  float theta12; // angle between the 2 gamma
  float energy;
  int ind_i;   // index in RecoParts of gamma i
  int ind_j;   // index in RecoParts of gamma j
  float delta; // difference with input mass
};

//===================================
//===================================
RVec<pi0_candidate>
build_pi0(const RVec<int> &ids,
          const RVec<edm4hep::ReconstructedParticleData> &in);

//============================
// return TLorentzVector of missing energy (can use custom RecoParts dataset VS
// MissingET by default needs all)
TLorentzVector missingTLV(float ecm,
                          const RVec<edm4hep::ReconstructedParticleData> &parts,
                          float p_cutoff = 0.0);

//===================================
// return RVec 2 of taus energy, solving for neutrinos in collinear
// approximation kinematics given by appendix in: "Measurement of Z' couplings
// at future hadron colliders through decays to r leptons"
RVec<float>
collinear_approx(const RVec<edm4hep::ReconstructedParticleData> &rps1,
                 const RVec<edm4hep::ReconstructedParticleData> &rps2,
                 double m_parent);

//===================================
// return RVec<MCParticleData> of MC particles from ids list
RVec<edm4hep::MCParticleData> getMC(const RVec<int> &ids,
                                    const RVec<edm4hep::MCParticleData> &in);

struct RPTruthInfo {
  TLorentzVector p4; // TLV PxPyPzM
  float rp_energy;   // particle energy
  int rp_index;      // Reco
  int mc_index;      // associated MCid
  int mc_pdg;
  float mc_energy;
  RVec<int> mc_parent_index;     // MC parent id
  RVec<int> mc_parent_pdg;       // MC parent pdg
  RVec<int> mc_parent_genStatus; // generator status of MC parent ()
  float mc_costheta_star;        // CosThetaStar for reweighting
  float costheta_thrust;         // costheta with respect to thrust
  int hemisphere; // hemisphere sign 0,1 NOT +1 -1!!! makes for cycles easy
  int charge;     // reco charge
  int mc_event;
  int reco_event;

  float mc_weight_plus;
  float mc_weight_minus;

  int ntot;
};

RVec<RPTruthInfo> buildRPTruthCollection(
    const RVec<edm4hep::ReconstructedParticleData> &reco,
    const RVec<edm4hep::MCParticleData> &mc,
    const RVec<int> &parents,     // Particle0
    const RVec<int> &rp2mc_index, // RP2MC_index
    const RVec<float> &costheta,
    const RVec<int> &sel_ids, // specific particle type ids
    const RVec<int> &MC_event,
    const RVec<int> &reco_event); // specific particle type ids

// generic functions
// ENERGY
RVec<float> truth_e(const RVec<RPTruthInfo> &truth);
// COSTHETA
RVec<float> truth_costheta(const RVec<RPTruthInfo> &truth);
// COSTHETA*
RVec<float> truth_costhetastar(const RVec<RPTruthInfo> &truth);
// PDG
RVec<int> truth_pdg(const RVec<RPTruthInfo> &truth);
// weight
RVec<float> truth_weight(const RVec<RPTruthInfo> &truth, bool foo);

// return vector of PDG int for given type of event
// truth sono le struct dei fotoni
RVec<int> parentPDG(const RVec<RPTruthInfo> &truth, const RVec<event> &ev,
                    const RVec<int> &MC_ev, int Nph, int type);

// return vector of energy for given type of event
RVec<float> pdgtype_e(const RVec<RPTruthInfo> &truth, const RVec<event> &ev,
                      const RVec<int> &MC_ev, int Nph, int type, bool foo_type,
                      int pdg, bool foo_pdg);

// return vector of energy for given type of event
RVec<float> confusion_e(const RVec<RPTruthInfo> &truth, int mc_event,
                        bool foomc, int reco_event, bool fooreco, int pdg,
                        bool foopdg, double P_cutoff);

RVec<float> true_e(const RVec<RPTruthInfo> &truth, int mc_event, bool foomc,
                   int reco_event, bool fooreco, int pdg, bool foopdg,
                   double P_cutoff);

// return vector of costheta (with respect to thrust) for given type of event
RVec<float> confusion_theta(const RVec<RPTruthInfo> &truth, int mc_event,
                            bool foomc, int reco_event, bool fooreco, int pdg,
                            bool foopdg);

// return vector of P for given type of event
RVec<float> confusion_p(const RVec<RPTruthInfo> &truth, int mc_event,
                        bool foomc, int reco_event, bool fooreco, int pdg,
                        bool foopdg, double P_cutoff);

// return energy of given mc particle (no .energy method is present)
float get_mc_e(const edm4hep::MCParticleData &mc);

// Function to calculate cos(theta*) using Truth MC 4-vectors
float GetCosThetaStar(const TLorentzVector &p4_tau_lab,
                      const TLorentzVector &p4_pi_lab);

// MC CLASSIFICATION
RVec<int> classify_mc_event(const RVec<edm4hep::MCParticleData> &mc,
                            const RVec<int> &daughters);
// auxiliary: pass tau id, MCParticle, Particle1 (daughters)
int mc_classification(int &tau, const RVec<edm4hep::MCParticleData> &mc,
                      const RVec<int> &ind);
// return pdg of tau daughters
RVec<int> print_dgt_pdg(const RVec<edm4hep::MCParticleData> &mc,
                        const RVec<int> &ind);
// return genstatus of [0] in get_decay
int print_tau_genstatus(const RVec<int> &tau,
                        const RVec<edm4hep::MCParticleData> &mc);

// return vector of histograms weights for given type of event
RVec<float> confusion_weight(const RVec<RPTruthInfo> &truth, int mc_event,
                             bool foomc, int reco_event, bool fooreco, int pdg,
                             bool foopdg, bool fooweight);

} // namespace Ztautau

#endif
