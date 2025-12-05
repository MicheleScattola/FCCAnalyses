// -*- C++ -*-
//
/** FCCAnalysis module: Z -> tau tau events
 *
 * \file Ztautau.cc
 * \author Michele Scattola <michele.scattola@studenti.unimi.it>
 */

#include "FCCAnalyses/Ztautau.h"

using namespace std;

namespace Ztautau {

//===================================
// THRUST
//===================================
// custom getThrustPointing using charge instead of energy
RVec<float> getThrustPointing(const RVec<float> &charge,
                              const RVec<float> &thrust,
                              const RVec<float> &costheta) {

  // copy
  RVec<float> out = thrust;

  // build total charges
  float c1 = 0.;
  float c2 = 0.;
  for (size_t i = 0; i < charge.size(); i++) {
    if (costheta[i] >= 0) {
      c1 += charge[i];
    } else if (costheta[i] < 0) {
      c2 += charge[i];
    }
  }

  // check c_pos and c_neg values. Adjust thrust direction accordingly
  // thrust starts in c1 direction
  bool flip = false;
  if (c1 > 0.0 && c2 <= 0.0) {
    flip = false;
  } else if (c1 <= 0.0 && c2 > 0.0) {
    flip = true;
  } else {
    // limit cases: both < or > 0 , flip towards "less negative" or "more
    // positive"
    if (c2 > c1)
      flip = true;
  }

  if (flip) {
    out[1] = -out[1];
    out[3] = -out[3];
    out[5] = -out[5];
  }

  return out;
}

//===================================
// SELECT PIONS
//===================================
// select pion indices as hadron with mass = 0.139570
RVec<int> sel_pions_id(const RVec<edm4hep::ReconstructedParticleData> &in,
                       const int charge) {
  RVec<int> out;

  for (int i = 0; i < in.size(); i++) {
    if (in[i].mass > 0.138 && in[i].mass < 0.14 && in[i].type == 0) {
      if (abs(in[i].charge) == charge) {
        out.push_back(i);
      }
    }
  }
  return out;
}

//===================================
// myEVENT
//===================================
// return event struct
RVec<myEvent> myget_event(const RVec<int> &mu_ids, const RVec<int> &el_ids,
                          const RVec<int> &pi_ids, const RVec<int> &ph_ids,
                          const RVec<edm4hep::ReconstructedParticleData> &rps,
                          const RVec<float> &rps_costheta, const bool masscheck,
                          const RVec<edm4hep::MCParticleData> &mc,
                          const RVec<int> &daughters) {

  // collect particles
  RVec<edm4hep::ReconstructedParticleData> mu_tot =
      ReconstructedParticle::get(mu_ids, rps);
  RVec<edm4hep::ReconstructedParticleData> el_tot =
      ReconstructedParticle::get(el_ids, rps);
  RVec<edm4hep::ReconstructedParticleData> pi_tot =
      ReconstructedParticle::get(pi_ids, rps);
  RVec<edm4hep::ReconstructedParticleData> ph_tot =
      ReconstructedParticle::get(ph_ids, rps);
  // collect their costheta vectors
  RVec<float> mu_costheta = get_elements_by_index(rps_costheta, mu_ids);
  RVec<float> el_costheta = get_elements_by_index(rps_costheta, el_ids);
  RVec<float> pi_costheta = get_elements_by_index(rps_costheta, pi_ids);
  RVec<float> ph_costheta = get_elements_by_index(rps_costheta, ph_ids);

  // cycle on hemisperhes, 0 is positive charge, 1 negative
  RVec<myEvent> out;
  bool hemisphere = true;

  for (int i = 0; i < 2; i++) {
    myEvent ev;

    // select particles in hemisphere
    // bool hemisphere starts as true for positive hemi, changing after first
    // loop to false for negative hemi
    RVec<edm4hep::ReconstructedParticleData> mu =
        ReconstructedParticle::sel_axis(hemisphere)(mu_costheta, mu_tot);
    RVec<edm4hep::ReconstructedParticleData> el =
        ReconstructedParticle::sel_axis(hemisphere)(el_costheta, el_tot);
    RVec<edm4hep::ReconstructedParticleData> pi =
        ReconstructedParticle::sel_axis(hemisphere)(pi_costheta, pi_tot);
    RVec<edm4hep::ReconstructedParticleData> ph =
        ReconstructedParticle::sel_axis(hemisphere)(ph_costheta, ph_tot);
    // fill struct
    ev.n_mu = mu.size();
    ev.n_el = el.size();
    ev.n_pi = pi.size();
    ev.n_ph = ph.size();

    // lambda helper to fill P4 and add energy & charge
    auto fill_collection = [&](const auto &input_particles,
                               RVec<TLorentzVector> &out_p4) {
      out_p4.reserve(input_particles.size());

      for (const auto &p : input_particles) {
        ev.m_RecoCharge += p.charge;
        ev.m_RecoEnergy += p.energy;

        TLorentzVector tlv;
        tlv.SetPxPyPzE(p.momentum.x, p.momentum.y, p.momentum.z, p.energy);
        out_p4.push_back(tlv);
      }
    };
    // fill collections
    fill_collection(mu, ev.m_muP4);
    fill_collection(el, ev.m_elP4);
    fill_collection(pi, ev.m_piP4);
    fill_collection(ph, ev.m_phP4);

    // RECO EVT CLASSIFICATION

    // Leptonic
    if ((ev.n_mu == 1 || ev.n_el == 1) && ev.n_pi == 0) {
      ev.m_type = classify_lep(ev);
      // TODO: ADD WEIGHT CALCULATION
      // TODO: ADD OPTIMAL VARIABLE CALCULATION
    }
    // Hadronic
    else if (ev.n_pi == 1 && ev.n_mu == 0 && ev.n_el == 0) {

      if (ev.m_piP4.size() > 0) {
        ev.m_pi_e = ev.m_piP4[0].E();
        ev.m_type = classify_pion(ev, masscheck);
        // TODO: ADD OPTIMAL VARIABLE CALCULATION
      } else {
        // debug
        cerr << "CRITICAL ERROR: n_pi is 1 but vector is empty inside loop!"
             << endl;
        ev.m_type = 0;
      }
    }
    // 3 prong
    else if (ev.n_pi == 3 && ev.n_mu == 0 && ev.n_el == 0 && ev.n_ph == 0) {
      // mass limit
      TLorentzVector p3pi = ev.m_piP4[0] + ev.m_piP4[1] + ev.m_piP4[2];
      if (masscheck && p3pi.M() < SM_TAU) {
        ev.m_type = 5; // Type 5: a1 (3-prong mode)
                       // TODO: ADD WEIGHT CALCULATION
                       // TODO: ADD OPTIMAL VARIABLE CALCULATION
      } else {
        ev.m_type = 0;
      }
    }

    // MC EVENT CLASSIFICATION & WEIGHTING

    // find tau, loop on daughters for event type, then weight based on evt type
    for (int i = 0; i < mc.size(); i++) {
      const auto &p = mc[i];
      // skip if not matching tau
      if (abs(p.PDG) != 15 || p.charge * ev.m_RecoCharge < 0 ||
          p.generatorStatus != 2)
        continue;
      // found tau with matching charge and decayed status
      ev.m_tauMCindex = i;
      TLorentzVector p4_tau_lab;
      p4_tau_lab.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
      // cycle daugthers and find event type
      int pb = p.daughters_begin;
      int pe = p.daughters_end;

      // sanity check
      if (pe == pb) {
        cerr << "[ERROR]: tau has no daughters in MC classification!" << endl;
        break;
      }
      // loop daughters, collect pdgs, classify
      RVec<int> dau_pdgs;
      for (int i = pb; i < pe; i++) {
        int dau_idx = daughters[i];
        const auto &dau = mc[dau_idx];
        dau_pdgs.push_back(abs(dau.PDG));
      }
      // classify
      ev.m_MCtype = 0;
      ev.m_MCtype = classify_MC(dau_pdgs);

      // weight calculation
      if (ev.m_MCtype == 3) {
        pion_weight(ev, mc, daughters);
      } else if (ev.m_MCtype == 4) {
        rho_weight(ev, mc, daughters);
      } else if (ev.m_MCtype == 5) {
        a1_weight(ev, mc, daughters);
      }

      // check on negative weights
      if (ev.m_MCweight_plus < 0.0 || ev.m_MCweight_minus < 0.0) {
        cerr << "[WARNING]: negative MC weight!" << endl;
      }
      // now exit the loop
      break;
    }

    // push back and change hemisphere
    out.push_back(ev);
    hemisphere = false;
  }

  return out;
}

// ==========================================
int classify_MC(const RVec<int> &pdgs) {
  int n_mu = 0, n_el = 0, n_pi = 0, n_pi0 = 0, n_ph = 0;
  for (const auto &p : pdgs) {
    if (p == 13)
      n_mu++;
    else if (p == 11)
      n_el++;
    // else if (p == 211 || p == 321 || p == 323) n_pi++; // do not distinguish
    // between pi or Kaon
    else if (p == 211)
      n_pi++; // do not distinguish between pi or Kaon
    else if (p == 111)
      n_pi0++;
    else if (p == 22)
      n_ph++;
  }
  // classification
  if (n_mu == 1 && n_el == 0 && n_pi == 0)
    return 1; // mu
  else if (n_mu == 0 && n_el == 1 && n_pi == 0)
    return 2; // el
  else if (n_mu == 0 && n_el == 0 && n_pi == 1 && n_pi0 == 0)
    return 3; // pi
  else if (n_mu == 0 && n_el == 0 && n_pi == 1 && n_pi0 == 1)
    return 4; // rho
  else if (n_mu == 0 && n_el == 0 && n_pi == 1 && n_pi0 == 2)
    return 5; // a1 (1prong)
  else if (n_mu == 0 && n_el == 0 && n_pi == 3)
    return 5; // a1 (3prong)

  return 0;
}

// ==========================================
int classify_lep(const myEvent &ev) {
  // Check MUON: 1 mu, 0 others
  if (ev.n_mu == 1 && ev.n_el == 0 && ev.n_pi == 0 && ev.n_ph == 0) {
    return 1; // Type 1: Muon
  }
  // Check ELECTRON: 1 el, 0 others
  if (ev.n_mu == 0 && ev.n_el == 1 && ev.n_pi == 0 && ev.n_ph == 0) {
    return 2; // Type 2: Electron
  }

  return 0;
}

// ==========================================
int classify_pion(const myEvent &ev, bool masscheck) {

  // assuming 1 pi and 0 leptons
  if (ev.m_piP4.empty()) {
    cerr << "ERROR in classify_pion: n_pi=" << ev.n_pi
         << " but vector is empty!" << endl;
    return 0;
  }
  // total visible P4
  TLorentzVector p4_vis = ev.m_piP4[0];
  for (const auto &ph_p4 : ev.m_phP4) {
    p4_vis += ph_p4;
  }
  // note: .M() returns invariant mass P4
  // also += method always returns a PxPyPzE vector (see ROOT docs)
  float mass_vis = p4_vis.M();

  // basic limit on tau mass
  if (masscheck && mass_vis > SM_TAU) {
    return 0;
  }

  // classification
  if (ev.n_ph == 0) {
    return 3; // Type 3: Single Pion
  } else if (ev.n_ph >= 1 && ev.n_ph <= 2) {
    return 4; // Type 4: Rho (pi + 1-2 gamma)
  } else if (ev.n_ph >= 3) {
    return 5; // Type 5: a1 -> pi + 2pi0 -> pi + 4gamma
  }

  return 0;
}

RVec<int> get_type_safe(const RVec<myEvent> &evs) {
  RVec<int> out;
  out.reserve(evs.size()); // Riserva memoria
  for (const auto &e : evs) {
    out.push_back(e.m_type);
  }
  return out;
}

RVec<float> get_energy_safe(const RVec<myEvent> &evs) {
  RVec<float> out;
  out.reserve(evs.size());
  for (const auto &e : evs) {
    // Usa il nome corretto della variabile nella struct
    out.push_back(e.m_pi_e);
  }
  return out;
}

// ==========================================
// RE-WEIGHTING FUNCTIONS
// ==========================================

float calc_Ptau(const TLorentzVector &p4_tau) {
  float costheta = p4_tau.CosTheta();
  float Ptau = -(SM_Atau * (1 + costheta * costheta) + 2 * SM_Atau * costheta) /
               (1 + costheta * costheta + 2 * SM_Atau * SM_Atau * costheta);
  return Ptau;
}

// Function to calculate cos(theta*) using Truth MC 4-vectors
float GetCosThetaStar(const TLorentzVector &p4_tau_lab,
                      const TLorentzVector &p4_pi_lab) {

  // BoostVector() returns beta of the particle (Lab -> Particle)
  TVector3 boost_to_rest = -p4_tau_lab.BoostVector();
  // copying and boosting
  TLorentzVector p4_pi_rest = p4_pi_lab;
  p4_pi_rest.Boost(boost_to_rest);

  // get tau's direction
  TVector3 tau_dir_lab = p4_tau_lab.Vect().Unit();

  // calculate angle between tau's flight and pi boosted in tau rest frame
  float angle = p4_pi_rest.Vect().Angle(tau_dir_lab);

  return cos(angle);
}

void pion_weight(myEvent &ev, const RVec<edm4hep::MCParticleData> &mc,
                 const RVec<int> &daughters) {

  const float charge = ev.m_RecoCharge;
  const int tau_idx = ev.m_tauMCindex;
  float z = 0.; // costheta star
  float alpha = 1.;
  float Ptau = 0.; // recalculated from tau p4

  // use tau index to find tau directly
  if (tau_idx < 0 || tau_idx >= mc.size()) {
    cerr << "[ERROR]: Invalid tau index" << endl;
    return;
  }
  // tau found
  const auto &p = mc[tau_idx];
  TLorentzVector p4_tau_lab;
  p4_tau_lab.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
  // cycle daugthers and find pion
  int pb = p.daughters_begin;
  int pe = p.daughters_end;

  Ptau = calc_Ptau(p4_tau_lab);

  for (int i = pb; i < pe; i++) {
    int dau_idx = daughters[i];
    const auto &dau = mc[dau_idx];
    // if (abs(dau.PDG) == 211 || abs(dau.PDG) == 321 || abs(dau.PDG) == 323) {
    if (abs(dau.PDG) == 211) {
      // found pion daughter
      ev.m_found = true;
      TLorentzVector p4_pi_lab;
      p4_pi_lab.SetXYZM(dau.momentum.x, dau.momentum.y, dau.momentum.z,
                        dau.mass);
      // boosting pion into tau rest frame
      z = GetCosThetaStar(p4_tau_lab, p4_pi_lab);
      // now exit the loops
      break;
    }
  }
  // weight
  float w_plus = (1 + alpha * z) / (1 + alpha * Ptau * z);
  float w_minus = (1 - alpha * z) / (1 - alpha * Ptau * z);

  ev.m_MCweight_plus = w_plus;
  ev.m_MCweight_minus = w_minus;
}

void rho_weight(myEvent &ev, const RVec<edm4hep::MCParticleData> &mc,
                const RVec<int> &daughters) {

  const float charge = ev.m_RecoCharge;
  const int tau_idx = ev.m_tauMCindex;
  float z = 0.;       // costheta star
  float alpha = 0.46; // recalculate
  float Ptau = 0.;

  // use tau index to find tau directly
  if (tau_idx < 0 || tau_idx >= mc.size()) {
    cerr << "[ERROR]: Invalid tau index" << endl;
    return;
  }
  // tau found
  const auto &p = mc[tau_idx];
  TLorentzVector p4_tau_lab;
  p4_tau_lab.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
  // cycle daugthers
  int pb = p.daughters_begin;
  int pe = p.daughters_end;

  Ptau = calc_Ptau(p4_tau_lab);
  TLorentzVector p4_rho_lab;

  for (int i = pb; i < pe; i++) {
    // loop on daughters, find resonance, store p4
    int dau_idx = daughters[i];
    const auto &dau = mc[dau_idx];
    // save pion
    if (abs(dau.PDG) == 211) {
      ev.m_found = true;
      TLorentzVector p4_pi_lab;
      p4_pi_lab.SetXYZM(dau.momentum.x, dau.momentum.y, dau.momentum.z,
                        dau.mass);
      p4_rho_lab += p4_pi_lab;
    } else if (abs(dau.PDG) == 111) {
      // save pi0
      TLorentzVector p4_pi0_lab;
      p4_pi0_lab.SetXYZM(dau.momentum.x, dau.momentum.y, dau.momentum.z,
                         dau.mass);
      p4_rho_lab += p4_pi0_lab;
    }
  }
  float mRho = p4_rho_lab.M();
  z = GetCosThetaStar(p4_tau_lab, p4_rho_lab);
  alpha =
      (SM_TAU * SM_TAU - 2 * mRho * mRho) / (SM_TAU * SM_TAU + 2 * mRho * mRho);
  // weights
  ev.m_MCweight_plus = (1 + alpha * z) / (1 + alpha * Ptau * z);
  ev.m_MCweight_minus = (1 - alpha * z) / (1 - alpha * Ptau * z);
}

void a1_weight(myEvent &ev, const RVec<edm4hep::MCParticleData> &mc,
               const RVec<int> &daughters) {

  const float charge = ev.m_RecoCharge;
  const int tau_idx = ev.m_tauMCindex;
  float z = 0.;       // costheta star
  float alpha = 0.12; // recalculate ?? worth it??
  float Ptau = 0.;

  // use tau index to find tau directly
  if (tau_idx < 0 || tau_idx >= mc.size()) {
    cerr << "[ERROR]: Invalid tau index" << endl;
    return;
  }
  // tau found
  const auto &p = mc[tau_idx];
  TLorentzVector p4_tau_lab;
  p4_tau_lab.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
  // cycle daugthers
  int pb = p.daughters_begin;
  int pe = p.daughters_end;

  Ptau = calc_Ptau(p4_tau_lab);
  TLorentzVector p4_rho_lab;

  for (int i = pb; i < pe; i++) {
    // loop on daughters, find resonance, store p4
    int dau_idx = daughters[i];
    const auto &dau = mc[dau_idx];
    // save pion
    if (abs(dau.PDG) == 211) {
      ev.m_found = true;
      TLorentzVector p4_pi_lab;
      p4_pi_lab.SetXYZM(dau.momentum.x, dau.momentum.y, dau.momentum.z,
                        dau.mass);
      p4_rho_lab += p4_pi_lab;
    } else if (abs(dau.PDG) == 111) {
      // save pi0
      TLorentzVector p4_pi0_lab;
      p4_pi0_lab.SetXYZM(dau.momentum.x, dau.momentum.y, dau.momentum.z,
                         dau.mass);
      p4_rho_lab += p4_pi0_lab;
    }
  }

  z = GetCosThetaStar(p4_tau_lab, p4_rho_lab);
  // weights
  ev.m_MCweight_plus = (1 + alpha * z) / (1 + alpha * Ptau * z);
  ev.m_MCweight_minus = (1 - alpha * z) / (1 - alpha * Ptau * z);
}

// ==========================================
// EXTRACT OPTIMAL VARIABLES
// ==========================================
RVec<float> get_lepton_e(const RVec<myEvent> &evs,
                                const int mc_type, const bool bool_mc,
                                const int reco_type, const bool bool_reco){

  RVec<float> out;                                 
  for(const auto &e : evs){
    // check mc event
    if(bool_mc && e.m_MCtype != mc_type) continue;
    // check reco event
    if(bool_reco && e.m_type != reco_type) continue;
    // leptonic x variable
    if(e.n_mu > 0) out.push_back(e.m_muP4[0].E());
    else if(e.n_el > 0) out.push_back(e.m_elP4[0].E());

  }
  return out;
};
// ==========================================
RVec<float> get_hadron_e(const RVec<myEvent> &evs,
                                const int mc_type, const bool bool_mc,
                                const int reco_type, const bool bool_reco){
  
  RVec<float> out;                                 
  for(const auto &e : evs){
    // check mc event
    if(bool_mc && e.m_MCtype != mc_type) continue;
    // check reco event
    if(bool_reco && e.m_type != reco_type) continue;
    // pion x variable
    if(n_pi > 0) out.push_back(e.m_piP4[0].E() / E_TAU);

  }
  return out;
};

RVec<int> get_pi_mask(const RVec<myEvent> &evs) {
    RVec<int> mask;
    for (const auto &e : evs) {
        // apply filter to pion events
        if (e.m_type == 3) {
            // safety check
            if (e.m_piP4.empty()) {
                mask.push_back(0); 
            } else {
                // energy cut
                if (e.m_piP4[0].E() > 2.0) mask.push_back(1); 
                else mask.push_back(0); 
            }
        } 
        // pass all other events
        else {
            mask.push_back(1);
        }
    }
    return mask;
}

/*
//===================================
// auxiliary function for classification
event classify(int n_mu, int n_el, int n_pi, int n_ph, const RVec<float> &ph_e,
               float ph_e_cutoff) {

  // collect max energy in photons collection
  event out;
  out.m_ph_e_max = 0.;
  if (ph_e.size() > 0) {
    out.m_ph_e_max = ROOT::VecOps::Max(ph_e);
    for (int i = 0; i < ph_e.size(); i++) {
      out.m_ph_e.push_back(ph_e[i]);
      out.m_sum_ph_e += ph_e[i];
    }
  }

  // assign struct values
  out.m_muN = n_mu;
  out.m_elN = n_el;
  out.m_piN = n_pi;
  out.m_phN = n_ph;

  // default: undefined
  out.type = 0;

  // 1: tau -> mu + else
  if (n_mu == 1 && n_el == 0 && n_pi == 0) {
    out.type = 1;
  }
  // 2: tau -> el + else
  else if (n_mu == 0 && n_el == 1 && n_pi == 0) {
    out.type = 2;
  }
  // 3: tau -> pi (no photons)
  else if (n_mu == 0 && n_el == 0 && n_pi == 1 && n_ph == 0) {
    out.type = 3;
  }

  // 4: tau -> rho (1 pi + 1-2 photons)
  else if (n_mu == 0 && n_el == 0 && n_pi == 1 && n_ph > 0 && n_ph <= 2) {
    out.type = 4;
  }
  // 5: tau -> a1 (3pi)
  else if (n_mu == 0 && n_el == 0 && n_pi == 3 && n_ph == 0) {
    out.type = 5;
  }
  // 5: tau -> a1 (1pi + >=3 gamma) , assumption we lose max 1ph
  else if (n_mu == 0 && n_el == 0 && n_pi == 1 && n_ph >= 3) {
    out.type = 5;
  }

  return out;

  // TODO: handle a1 cases
}

//===================================
// return event struct
RVec<event> get_event(const RVec<edm4hep::ReconstructedParticleData> &mu,
                      const RVec<float> &mu_costheta,
                      const RVec<edm4hep::ReconstructedParticleData> &el,
                      const RVec<float> &el_costheta,
                      const RVec<edm4hep::ReconstructedParticleData> &pi,
                      const RVec<float> &pi_costheta,
                      const RVec<edm4hep::ReconstructedParticleData> &ph,
                      const RVec<float> &ph_costheta, const float ph_e_cutoff) {

  RVec<int> mu_n = aux_n(mu_costheta);
  RVec<int> el_n = aux_n(el_costheta);
  RVec<int> ph_n = aux_n(ph_costheta);
  RVec<int> pi_n = aux_n(pi_costheta);

  RVec<float> ph_pos_e = aux_e(0, ph, ph_costheta);
  RVec<float> ph_neg_e = aux_e(1, ph, ph_costheta);

  event classified_pos =
      classify(mu_n[0], el_n[0], pi_n[0], ph_n[0], ph_pos_e, ph_e_cutoff);
  event classified_neg =
      classify(mu_n[1], el_n[1], pi_n[1], ph_n[1], ph_neg_e, ph_e_cutoff);

  RVec<event> out{classified_pos, classified_neg};

  return out;
}

//===================================
// auxiliary function, return RVec of particles energies in passed hemisphere
// sign useful for photons
RVec<int> aux_n(const RVec<float> &legs_costheta) {

  int pos = 0;
  int neg = 0;
  for (int i = 0; i < legs_costheta.size(); i++) {
    if (legs_costheta[i] >= 0)
      pos++;
    else
      neg++;
  }

  return RVec<int>{pos, neg};
}

//===================================
// auxiliary function, return RVec of particles energies in passed hemisphere
// sign
RVec<float> aux_e(int hemisphere,
                  const RVec<edm4hep::ReconstructedParticleData> &legs,
                  const RVec<float> &legs_costheta) {

  // auxiliary bool for hemi sign
  bool hemi;
  if (hemisphere == 0)
    hemi = true;
  if (hemisphere == 1)
    hemi = false;
  // retrieve corresponding reconstructed particles
  RVec<edm4hep::ReconstructedParticleData> sel;
  sel = ReconstructedParticle::sel_axis(hemi)(legs_costheta, legs);

  RVec<float> out = ReconstructedParticle::get_e(sel);

  return out;
}

//===================================
// montecarlo classification
int MC_classified(const RVec<int> &TauMu, const RVec<int> &TauEl,
                  const RVec<int> &TauPi, const RVec<int> &TauRho) {
  if (TauMu.size() > 0)
    return 1;
  else if (TauEl.size() > 0)
    return 2;
  else if (TauPi.size() > 0)
    return 3;
  else if (TauRho.size() > 0)
    return 4;
  else
    return 0;
}

// vector of photon energies given conditions
// bool for TrueMC rho event
// RVec<event> &ev are RecoPart properties of event (per hemisphere)
// int N number of RecoPhotons collected
// in type is event type 0,1,2,3,4 ..
RVec<float> study_ph(const RVec<int> &MC_event, const RVec<event> &ev, bool foo,
                     int N_ph, bool wrong_events, int type) {

  RVec<float> e;

  // auxiliary lambda to fill input vector with event's ph_e
  auto aux = [](const event P, RVec<float> &in) {
    for (int i = 0; i < P.m_ph_e.size(); i++)
      in.push_back(P.m_ph_e[i]);
  };
  // cycle on hemispheres
  for (int i = 0; i < 2; i++) {
    // if looking for TrueMC rho
    if (foo && !wrong_events) {
      if (MC_event[i] == type && ev[i].m_piN == 1 && ev[i].m_phN == N_ph)
        aux(ev[i], e);
    } else if (!foo && !wrong_events) {
      if (ev[i].m_piN == 1 && ev[i].m_phN == N_ph)
        aux(ev[i], e);
    } else if (wrong_events) {
      // looking for NON-RHO EVENTS which get reconstructed with 1pi+N_ph
      // photons (mainly rho&a1)
      if (MC_event[i] != type && ev[i].m_piN == 1 && ev[i].m_phN == N_ph)
        aux(ev[i], e);
    }
  }
  return e;
}

// return sum of energies
RVec<float> study_ph_sum(const RVec<int> &MC_event, const RVec<event> &ev,
                         bool foo, int N_ph, bool wrong_events, int type) {

  RVec<float> e; // needed if both hemispheres have same event

  // auxiliary lambda to fill input vector with event's ph_e
  // cycle on hemispheres
  for (int i = 0; i < 2; i++) {
    // if looking for TrueMC rho
    if (foo && !wrong_events) {
      if (MC_event[i] == type && ev[i].m_piN == 1 && ev[i].m_phN == N_ph)
        e.push_back(ev[i].m_sum_ph_e);
    } else if (!foo && !wrong_events) {
      if (ev[i].m_piN == 1 && ev[i].m_phN == N_ph)
        e.push_back(ev[i].m_sum_ph_e);
    } else if (wrong_events) {
      // looking for NON-PI0 EVENTS which get reconstructed with 1pi+N_ph
      // photons (mainly rho&a1)
      if (MC_event[i] != type && ev[i].m_piN == 1 && ev[i].m_phN == N_ph)
        e.push_back(ev[i].m_sum_ph_e);
    }
  }
  return e;
}
*/

/*

//===================================
// building all pairs for pi0 resonance from photons
RVec<edm4hep::ReconstructedParticleData>
pi0_resonance_pairs(const RVec<edm4hep::ReconstructedParticleData> &legs,
                    const float target_mass) {

  RVec<edm4hep::ReconstructedParticleData> out;
  const int n = legs.size();
  if (n < 2)
    return out;
  // loop on photons and construct resonant mass
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      TLorentzVector a, b, sum;
      a.SetPxPyPzE(legs[i].momentum.x, legs[i].momentum.y, legs[i].momentum.z,
                   legs[i].energy);
      b.SetPxPyPzE(legs[j].momentum.x, legs[j].momentum.y, legs[j].momentum.z,
                   legs[j].energy);
      sum = a + b;
      edm4hep::ReconstructedParticleData r{};
      r.charge = legs[i].charge + legs[j].charge;
      r.momentum.x = sum.Px();
      r.momentum.y = sum.Py();
      r.momentum.z = sum.Pz();
      r.mass = sum.M();
      r.energy = sum.E();
      out.emplace_back(r);
    }
  }

  // sort by difference to input mass
  sort(out.begin(), out.end(),
            [target_mass](const auto &A, const auto &B) {
              return fabs(target_mass - A.mass) <
                     fabs(target_mass - B.mass);
            });

  return out;
}

//===================================
RVec<pi0_candidate>
build_pi0(const RVec<int> &ids,
          const RVec<edm4hep::ReconstructedParticleData> &in) {

  RVec<pi0_candidate> out;
  int n = ids.size();
  // construct resonant mass from combinations of i,j photons
  for (int a = 0; a < n; a++) {
    for (int b = a + 1; b < n; b++) {

      // photons indices
      int i = ids[a];
      int j = ids[b];

      float Ei = in[i].energy;
      float Ej = in[j].energy;
      float Pxi = in[i].momentum.x;
      float Pyi = in[i].momentum.y;
      float Pzi = in[i].momentum.z;
      float Pxj = in[j].momentum.x;
      float Pyj = in[j].momentum.y;
      float Pzj = in[j].momentum.z;

      const float pi2 = Pxi * Pxi + Pyi * Pyi + Pzi * Pzi;
      const float pj2 = Pxj * Pxj + Pyj * Pyj + Pzj * Pzj;
      const float pi = sqrt(max(0.f, pi2));
      const float pj = sqrt(max(0.f, pj2));
      float dot = Pxi * Pxj + Pyi * Pyj + Pzi * Pzj;
      float costheta = dot / (pi * pj);

      pi0_candidate P;

      P.mass = sqrt(2. * Ei * Ej * (1. - (costheta)));
      P.theta12 = acos(costheta);
      P.energy = Ei + Ej;
      P.ind_i = i;
      P.ind_j = j;
      P.delta = abs(P.mass - P.mpi0);
      out.push_back(P);
    }
  }

  // sort by difference to pi0 mass
  sort(out.begin(), out.end(),
            [](const pi0_candidate &A, const pi0_candidate &B) {
              return A.delta < B.delta;
            });

  return out;
}

//===================================
// return TLorentzVector of missing energy from each event
TLorentzVector missingTLV(float ecm,
                          const RVec<edm4hep::ReconstructedParticleData> &parts,
                          float p_cutoff) {
  float px = 0, py = 0, pz = 0, e = 0;
  for (auto &p : parts) {
    if (sqrt(p.momentum.x * p.momentum.x + p.momentum.y * p.momentum.y) <
        p_cutoff)
      continue;
    px += -p.momentum.x;
    py += -p.momentum.y;
    pz += -p.momentum.z;
    e += p.energy;
  }
  TLorentzVector mis;
  mis.SetPxPyPzE(px, py, pz, ecm - e);

  return mis;
}

//===================================
// return RVec 2 of taus energy, solving for neutrinos in collinear
// approximation kinematics given by appendix in: "Measurement of Z' couplings
// at future hadron colliders through decays to r leptons"
RVec<float>
collinear_approx(const RVec<edm4hep::ReconstructedParticleData> &rps1,
                 const RVec<edm4hep::ReconstructedParticleData> &rps2,
                 double m_parent) {

  RVec<float> out;
  float p1x = 0., p1y = 0., p1z = 0., Evis1 = 0.;
  float p2x = 0., p2y = 0., p2z = 0., Evis2 = 0.;

  // collect variables in the 2 hemispheres
  for (int i = 0; i < rps1.size(); i++) {
    p1x += rps1[i].momentum.x;
    p1y += rps1[i].momentum.y;
    p1z += rps1[i].momentum.z;
    Evis1 += rps1[i].energy;
  }
  TLorentzVector tlv1;
  tlv1.SetPxPyPzE(p1x, p1y, p1z, Evis1);

  for (int i = 0; i < rps2.size(); i++) {
    p2x += rps2[i].momentum.x;
    p2y += rps2[i].momentum.y;
    p2z += rps2[i].momentum.z;
    Evis2 += rps2[i].energy;
  }
  TLorentzVector tlv2;
  tlv2.SetPxPyPzE(p2x, p2y, p2z, Evis2);

  // collinear approx
  float A = (2.0 * tlv1.Dot(tlv2)) / (m_parent * m_parent);
  float B = tlv1.Perp() / tlv2.Perp();
  float x1 = sqrt(A * B);
  float x2 = sqrt(A / B);

  out = {Evis1 / x1, Evis2 / x2};
  return out;
}

//===================================
// return RVec<MCParticleData> for pi0->gammagamma events
RVec<edm4hep::MCParticleData> getMC(const RVec<int> &ids,
                                    const RVec<edm4hep::MCParticleData> &in) {

  RVec<edm4hep::MCParticleData> out;

  RVec<int> idx = {ids[1], ids[2]};
  const int n = in.size();
  for (int i : idx) {
    if (i >= 0 && i < n && in[i].PDG == 22)
      out.push_back(in[i]);
  }

  return out;
}
*/
//===================================
//===================================
// CUSTOM DATA COLLECTION
//===================================
//===================================

// building custom struct for any reco particle type: mu,el,pi,gamma
RVec<RPTruthInfo> buildRPTruthCollection(
    const RVec<edm4hep::ReconstructedParticleData> &reco,
    const RVec<edm4hep::MCParticleData> &mc,
    const RVec<int> &parents,     // Particle0
    const RVec<int> &rp2mc_index, // RP2MC_index
    const RVec<float> &costheta,  // costheta with respect to thrust
    const RVec<int> &sel_ids,     // specific particle type ids
    const RVec<int> &MC_event, const RVec<int> &reco_event) {

  RVec<RPTruthInfo> out;
  out.reserve(sel_ids.size());

  const int nReco = reco.size();
  const int nMC = mc.size();

  for (int rp_idx : sel_ids) {
    if (rp_idx < 0 || rp_idx >= nReco)
      continue;

    RPTruthInfo info;

    // info from ReconstructedParticles
    const auto &rp = reco[rp_idx];
    TLorentzVector p4;
    p4.SetXYZM(rp.momentum.x, rp.momentum.y, rp.momentum.z, rp.mass);
    info.rp_energy = rp.energy;
    info.p4 = p4;
    info.rp_index = rp_idx;
    info.charge = rp.charge;
    info.costheta_thrust =
        (rp_idx < (int)costheta.size()) ? costheta[rp_idx] : -999.f;

    if (info.costheta_thrust >= 0.f) {
      info.hemisphere = 0;
      info.mc_event = MC_event[0];
      info.reco_event = reco_event[0];
    } else if (info.costheta_thrust < 0.f && info.costheta_thrust > -2.f) {
      info.hemisphere = 1;
      info.mc_event = MC_event[1];
      info.reco_event = reco_event[1];
    } else {
      info.hemisphere = -999;
      info.mc_event = -1;
      info.reco_event = -1;
    }

    // Matching MC information

    // exploiting  rp2mc map. Return -1 if bad case (it should already be -1
    // btw)
    int mc_idx = (rp_idx < (int)rp2mc_index.size()) ? rp2mc_index[rp_idx] : -1;
    info.mc_index = mc_idx;

    RVec<int> parent_idx;
    RVec<int> parent_pdg;
    RVec<int> parent_genStatus;
    if (mc_idx >= 0 && mc_idx < nMC) {
      const auto &pmc = mc[mc_idx];
      info.mc_pdg = pmc.PDG;
      info.mc_energy = get_mc_e(pmc);

      int pb = pmc.parents_begin;
      int pe = pmc.parents_end;
      // handle all non empty cases
      if (pb != pe) {
        // cycle on all possible parents (1 or more in strange cases?)
        for (int i = pb; i < pe; i++) {
          // corresponding id in the parents collection
          int temp = parents[i];
          parent_idx.push_back(temp);
          parent_pdg.push_back(mc[temp].PDG);
          parent_genStatus.push_back(mc[temp].generatorStatus);

          // collecting info for costheta star between meson and tau
          // check for pi and tau parent
          if (abs(mc[temp].PDG) == 15 &&
              (abs(pmc.PDG) == 211 || abs(pmc.PDG) == 321)) {
            TLorentzVector p4_tau_lab;
            p4_tau_lab.SetXYZM(mc[temp].momentum.x, mc[temp].momentum.y,
                               mc[temp].momentum.z, mc[temp].mass);
            TLorentzVector p4_pi_lab;
            p4_pi_lab.SetXYZM(pmc.momentum.x, pmc.momentum.y, pmc.momentum.z,
                              pmc.mass);

            // boosting pion into tau rest frame
            float z = GetCosThetaStar(p4_tau_lab, p4_pi_lab);
            info.mc_costheta_star = z;

          } else
            info.mc_costheta_star = -999.0f;
        }
      } else {
        // empty cases
        parent_idx = {-1};
        parent_pdg = {0};
      }
    }
    info.mc_parent_index = parent_idx;
    info.mc_parent_pdg = parent_pdg;
    info.mc_parent_genStatus = parent_genStatus;
    info.ntot = parent_idx.size();

    // reweighting
    float P = -0.150;
    float z = info.mc_costheta_star;
    float alpha = 1.;
    // coeff for rho
    if (info.mc_event == 4)
      alpha = 0.46;
    // the 'other' events are considered a1
    else if (info.mc_event == 5)
      alpha = 0.12;

    float B = alpha * z;
    if (B > 1.0)
      B = 1.0;
    if (B < -1.0)
      B = -1.0;

    info.mc_weight_plus = (1 + B) / (1 + B * P);
    info.mc_weight_minus = (1 - B) / (1 + B * P);
    if (z < -998.0f) {
      info.mc_weight_plus = 0;
      info.mc_weight_minus = 0;
    }

    out.push_back(info);
  }

  return out;
}

// generic functions
// ENERGY
RVec<float> truth_e(const RVec<RPTruthInfo> &truth) {

  RVec<float> out;
  for (const auto p : truth) {
    out.push_back(p.rp_energy);
  }

  return out;
}
// COSTHETA
RVec<float> truth_costheta(const RVec<RPTruthInfo> &truth) {

  RVec<float> out;
  for (const auto p : truth) {
    out.push_back(p.costheta_thrust);
  }

  return out;
}

// COSTHETA*
RVec<float> truth_costhetastar(const RVec<RPTruthInfo> &truth) {

  RVec<float> out;
  for (const auto p : truth) {
    out.push_back(p.mc_costheta_star);
  }

  return out;
}
// PDG
RVec<int> truth_pdg(const RVec<RPTruthInfo> &truth) {

  RVec<int> out;
  for (const auto p : truth) {
    out.push_back(p.mc_pdg);
  }

  return out;
}
// weight
RVec<float> truth_weight(const RVec<RPTruthInfo> &truth, bool foo) {

  RVec<float> out;
  for (const auto p : truth) {
    if (p.mc_costheta_star < -1)
      continue; // cases with -999.f
    if (foo)
      out.push_back(p.mc_weight_plus);
    else
      out.push_back(p.mc_weight_minus);
  }

  return out;
}

// return vector of PDG int for given type of event
// truth sono le struct dei fotoni
RVec<int> parentPDG(const RVec<RPTruthInfo> &truth, const RVec<event> &ev,
                    const RVec<int> &MC_ev, int Nph, int type) {

  RVec<int> out;
  RVec<int> unique_parents; // auxiliary for no duplicates

  for (int hemi = 0; hemi < 2; hemi++) {
    // selection
    if (MC_ev[hemi] == type)
      continue;
    if (ev[hemi].m_piN != 1)
      continue;
    if (ev[hemi].m_phN != Nph)
      continue;

    // cycle on truth
    for (const auto &ph : truth) {

      if (ph.hemisphere != hemi)
        continue;

      int n = ph.mc_parent_index.size();
      // store all pdg
      for (int j = 0; j < n; j++) {
        // no duplicates, check via same mc_parent_id
        bool already = false;
        int tempid = ph.mc_parent_index[j];
        for (const auto &u : unique_parents) {
          if (u == tempid) {
            already = true;
            break;
          }
        }
        // after check if its unique save data
        unique_parents.push_back(tempid);
        out.push_back(abs(ph.mc_parent_pdg[j]));
      }
    }
  }

  return out;
}

// return vector of ph_e for given type of event
// 2 boolens to allow custom request of matching event type and/or matching
// parent pdg
RVec<float> pdgtype_ph_e(const RVec<RPTruthInfo> &truth, const RVec<event> &ev,
                         const RVec<int> &MC_ev, int Nph, int type,
                         bool foo_type, int pdg, bool foo_pdg) {

  RVec<float> out;

  for (int hemi = 0; hemi < 2; hemi++) {
    // selection
    if (foo_type && MC_ev[hemi] != type)
      continue;
    // option to check
    if (ev[hemi].m_piN != 1)
      continue;
    if (ev[hemi].m_phN != Nph)
      continue;

    // cycle on truth
    for (const auto &ph : truth) {

      if (ph.hemisphere != hemi)
        continue;

      // check if ph is from ISR
      if (foo_pdg && abs(ph.mc_parent_pdg[0]) != pdg)
        continue;
      out.push_back(ph.rp_energy);
    }
  }

  return out;
}

// return vector of energy for given type of event
RVec<float> confusion_e(const RVec<RPTruthInfo> &truth, int mc_event,
                        bool foomc, int reco_event, bool fooreco, int pdg,
                        bool foopdg, double P_cutoff) {

  RVec<float> out;

  for (const auto p : truth) {
    if (foomc) {
      if (mc_event >= 0) {
        if (p.mc_event != mc_event)
          continue;
      } else {
        // looking for background events
        if (p.mc_event == -mc_event)
          continue;
      }
    }
    if (fooreco && p.reco_event != reco_event)
      continue;
    if (foopdg && abs(p.mc_pdg) != pdg)
      continue;
    if (p.p4.P() < P_cutoff)
      continue;

    out.push_back(p.rp_energy);
  }

  return out;
}

RVec<float> true_e(const RVec<RPTruthInfo> &truth, int mc_event, bool foomc,
                   int reco_event, bool fooreco, int pdg, bool foopdg,
                   double P_cutoff) {

  RVec<float> out;

  for (const auto p : truth) {
    if (foomc) {
      if (mc_event >= 0) {
        if (p.mc_event != mc_event)
          continue;
      } else {
        // looking for background events
        if (p.mc_event == -mc_event)
          continue;
      }
    }
    if (fooreco && p.reco_event != reco_event)
      continue;
    if (foopdg && abs(p.mc_pdg) != pdg)
      continue;
    if (p.p4.P() < P_cutoff)
      continue;

    out.push_back(p.mc_energy);
  }

  return out;
}

// return vector of costheta for given type of event
RVec<float> confusion_theta(const RVec<RPTruthInfo> &truth, int mc_event,
                            bool foomc, int reco_event, bool fooreco, int pdg,
                            bool foopdg) {

  RVec<float> out;

  for (const auto p : truth) {
    if (foomc) {
      if (mc_event >= 0) {
        if (p.mc_event != mc_event)
          continue;
      } else {
        // looking for background events
        if (p.mc_event == -mc_event)
          continue;
      }
    }
    if (fooreco && p.reco_event != reco_event)
      continue;
    if (foopdg && abs(p.mc_pdg) != pdg)
      continue;

    int sign;
    if (p.hemisphere == 0)
      sign = +1;
    else if (p.hemisphere == 1)
      sign = -1;
    // normalize angle with respect to thrust direction
    float norm_costheta = p.costheta_thrust * sign;
    out.push_back(fabs(norm_costheta));
  }

  return out;
}

// return vector of P for given type of event
RVec<float> confusion_p(const RVec<RPTruthInfo> &truth, int mc_event,
                        bool foomc, int reco_event, bool fooreco, int pdg,
                        bool foopdg, double P_cutoff) {

  RVec<float> out;

  for (const auto p : truth) {
    if (foomc) {
      if (mc_event >= 0) {
        if (p.mc_event != mc_event)
          continue;
      } else {
        // looking for background events
        if (p.mc_event == -mc_event)
          continue;
      }
    }
    if (fooreco && p.reco_event != reco_event)
      continue;
    if (foopdg && abs(p.mc_pdg) != pdg)
      continue;
    if (p.p4.P() < P_cutoff)
      continue;

    out.push_back(p.p4.P());
  }

  return out;
}

// return energy of given mc particle (no .energy method is present)
float get_mc_e(const edm4hep::MCParticleData &mc) {
  TLorentzVector tlv;
  tlv.SetXYZM(mc.momentum.x, mc.momentum.y, mc.momentum.z, mc.mass);

  return tlv.E();
}

//===================================
//===================================
// MC CLASSIFICATION
//===================================
//===================================

// Definition of Event Types
// 0: Other
// 1: Muon (tau -> mu nu nu)
// 2: Electron (tau -> e nu nu)
// 3: Pion (tau -> pi nu)
// 4: Rho (tau -> pi pi0 nu)
// 5: a1 (tau -> 3pi nu)

RVec<int> classify_mc_event(const RVec<edm4hep::MCParticleData> &mc,
                            const RVec<int> &daughters) {

  // results of type: {class_tau_plus, class_tau_minus}
  // default to 0
  RVec<int> classifications = {0, 0};

  // finding taus with GenStatus = 2

  int idx_tau_minus = -1;
  int idx_tau_plus = -1;

  for (size_t i = 0; i < mc.size(); ++i) {
    // pdg 15 and genstatus 2
    if (abs(mc[i].PDG) != 15)
      continue;
    if (mc[i].generatorStatus != 2)
      continue;

    if (mc[i].PDG == 15)
      idx_tau_minus = i;
    if (mc[i].PDG == -15)
      idx_tau_plus = i;
  }

  // Classify Tau Plus
  classifications[0] = mc_classification(idx_tau_plus, mc, daughters);

  // Classify Tau Minus
  classifications[1] = mc_classification(idx_tau_minus, mc, daughters);

  return classifications;
}

// auxiliary: pass tau id, MCParticle, Particle1 (daughters)
int mc_classification(int &tau, const RVec<edm4hep::MCParticleData> &mc,
                      const RVec<int> &ind) {

  int db = mc.at(tau).daughters_begin;
  int de = mc.at(tau).daughters_end;
  if (db == de)
    return -1; // particle is stable

  int n_mu = 0, n_el = 0, n_pi_ch = 0, n_pi_neu = 0, n_ph = 0;
  // loop on daughters
  for (int i = db; i < de; i++) {
    int temp = ind[i];
    int pdg = abs(mc[temp].PDG);

    if (pdg == 11)
      n_el++;
    else if (pdg == 13)
      n_mu++;
    else if (pdg == 211 || pdg == 321)
      n_pi_ch++; // kaons are grouped together
    else if (pdg == 111)
      n_pi_neu++;
    else if (pdg == 22)
      n_ph++;
  }

  // classification
  if (n_mu == 1)
    return 1;
  else if (n_el == 1)
    return 2;
  else if (n_pi_ch == 1 && n_pi_neu == 1)
    return 4;
  else if (n_pi_ch == 1 && n_pi_neu == 0)
    return 3;
  else if (n_pi_ch == 3)
    return 5;
  else if (n_pi_ch == 1 && n_pi_neu == 2)
    return 5;
  else
    return 0;
}

// return pdg of tau daughters
RVec<int> print_dgt_pdg(const RVec<edm4hep::MCParticleData> &mc,
                        const RVec<int> &ind) {
  RVec<int> out;
  int idx_tau_minus = -1;
  int idx_tau_plus = -1;

  for (size_t i = 0; i < mc.size(); ++i) {
    // pdg 15 and genstatus 2
    if (abs(mc[i].PDG) != 15)
      continue;
    if (mc[i].generatorStatus != 2)
      continue;

    if (mc[i].PDG == 15)
      idx_tau_minus = i;
    if (mc[i].PDG == -15)
      idx_tau_plus = i;
  }

  int db = mc.at(idx_tau_plus).daughters_begin;
  int de = mc.at(idx_tau_plus).daughters_end;

  // loop on daughters
  for (int i = db; i < de; i++) {
    int temp = ind[i];
    int pdg = mc[temp].PDG;
    out.push_back(pdg);
  }
  db = mc.at(idx_tau_minus).daughters_begin;
  de = mc.at(idx_tau_minus).daughters_end;
  // loop on daughters
  for (int i = db; i < de; i++) {
    int temp = ind[i];
    int pdg = mc[temp].PDG;
    out.push_back(pdg);
  }

  return out;
}

int print_tau_genstatus(const RVec<int> &tau,
                        const RVec<edm4hep::MCParticleData> &mc) {
  int temp = tau[0];
  return mc[temp].generatorStatus;
}

// return vector of histograms weights for given type of event
RVec<float> confusion_weight(const RVec<RPTruthInfo> &truth, int mc_event,
                             bool foomc, int reco_event, bool fooreco, int pdg,
                             bool foopdg, bool fooweight) {

  RVec<float> out;

  for (const auto p : truth) {
    if (foomc) {
      if (mc_event >= 0) {
        if (p.mc_event != mc_event)
          continue;
      } else {
        // looking for background events
        if (p.mc_event == -mc_event)
          continue;
      }
    }
    if (fooreco && p.reco_event != reco_event)
      continue;
    if (foopdg && abs(p.mc_pdg) != pdg)
      continue;

    if (fooweight)
      out.push_back(p.mc_weight_plus);
    if (!fooweight)
      out.push_back(p.mc_weight_minus);
  }

  return out;
}

} // namespace Ztautau
