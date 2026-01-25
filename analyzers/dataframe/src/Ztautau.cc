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
RVec<double> getThrustPointing(const RVec<double> &charge,
                              const RVec<double> &thrust,
                              const RVec<double> &costheta) {

  // copy
  RVec<double> out = thrust;

  // build total charges
  double c1 = 0.;
  double c2 = 0.;
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
                          const RVec<double> &rps_costheta,
                          const RVec<edm4hep::MCParticleData> &mc,
                          const RVec<int> &daughters,
                          const RVec<int> &rp2mc_idx,
                          const double &thrust_costheta,
                          const double &thrust_phi) {

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
  RVec<double> mu_costheta = get_elements_by_index(rps_costheta, mu_ids);
  RVec<double> el_costheta = get_elements_by_index(rps_costheta, el_ids);
  RVec<double> pi_costheta = get_elements_by_index(rps_costheta, pi_ids);
  RVec<double> ph_costheta = get_elements_by_index(rps_costheta, ph_ids);


  // cycle on hemisperhes, 0 is positive charge, 1 negative
  RVec<myEvent> out;
  bool hemisphere = true;

  for (int i = 0; i < 2; i++) {
    myEvent ev;
    ev.thrust_costheta = thrust_costheta;
    ev.thrust_phi = thrust_phi;
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
    // get ids vector for RP2MC
    RVec<int> mu_idx = get_idx(hemisphere,mu_costheta,mu_ids);
    RVec<int> el_idx = get_idx(hemisphere,el_costheta,el_ids);
    RVec<int> pi_idx = get_idx(hemisphere,pi_costheta,pi_ids);
    RVec<int> ph_idx = get_idx(hemisphere,ph_costheta,ph_ids);
    // fill struct
    ev.n_mu = mu.size();
    ev.n_el = el.size();
    ev.n_pi = pi.size();
    ev.n_ph = ph.size();

    // fill collections
    ev.m_RecoEnergy = 0.;
    ev.m_RecoMass = 0.;
    ev.m_RecoCharge = 0.;
    fill_collection(ev, mu, ev.m_muP4);
    fill_collection(ev, el, ev.m_elP4);
    fill_collection(ev, pi, ev.m_piP4);
    fill_collection(ev, ph, ev.m_phP4);

    // invariant mass check
    TLorentzVector p4_tot;
    for(auto &p : ev.m_muP4) p4_tot += p;
    for(auto &p : ev.m_elP4) p4_tot += p;
    for(auto &p : ev.m_piP4) p4_tot += p;
    for(auto &p : ev.m_phP4) p4_tot += p;
    
    // generic mass limit less than 2 GeV
    ev.m_RecoMass = p4_tot.M(); 
    if (ev.m_RecoMass > 2){
      ev.m_debug_mass = 1; // high mass
    }
    // mass limits on hadronic decays are done in classification

    // RECO EVT CLASSIFICATION
    // Leptonic
    if ( (ev.n_mu == 1 || ev.n_el == 1) && ev.n_pi == 0 ) {
      ev.m_type = classify_lep(ev,mu_idx,el_idx,rp2mc_idx);
    }
    // Hadronic
    else if (ev.n_pi > 0 && ev.n_mu == 0 && ev.n_el == 0) {
      ev.m_type = classify_pion(ev,pi_idx,rp2mc_idx);
    }
    
    else if ( (ev.n_mu == 1 || ev.n_el == 1) && ev.n_pi != 0 ) {
      ev.m_debug = 33;
    }

    // MC EVENT CLASSIFICATION & WEIGHTING

    // find tau, loop on daughters for event type, then weight based on evt type
    for (int i = 0; i < mc.size(); i++) {

      const auto &p = mc[i];
      // skip if not matching decayed tau
      // choose tau charge based on charge of hemisphere, in order they should be : POS , NEG
      // DO NOT RELY ON RECO CHARGES, BIASED!!!
      double temp_charge = 0.0;
      if(hemisphere) temp_charge = 1.0;
      else if (!hemisphere) temp_charge = -1.0;

      /*if (abs(p.PDG) != 15 || p.charge * temp_charge < 0 ||
          p.generatorStatus != 2)
        continue;*/
      if(abs(p.PDG) != 15 || p.charge * ev.m_RecoCharge < 0 ||
          p.generatorStatus != 2)
        continue;

  
      // cycle daugthers and find event type
      int pb = p.daughters_begin;
      int pe = p.daughters_end;
      bool tau_not_final = false;

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

        if(abs(dau.PDG) == 15) {
          tau_not_final = true;
          // exit
          break;
        }
      }
      // if tau has tau daughter skip this particle, it's not final tau
      if(tau_not_final) continue;

      // found tau with matching charge and decayed status
      ev.mc_tau_index = i;
      TLorentzVector p4_tau_lab;
      p4_tau_lab.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
      ev.mc_tauP4 = p4_tau_lab;

      ev.mc_daughters = dau_pdgs;
      // classify
      ev.mc_type = classify_MC(dau_pdgs);
      
      // weight calculation
      // weights also set trueMC daughter p4 and mass
      if (ev.mc_type == 1 || ev.mc_type == 2 ) {
      	lepton_weight(ev, mc, daughters);
      }	else if (ev.mc_type == 3) {
        pion_weight(ev, mc, daughters);
      } else if (ev.mc_type == 4) {
        rho_weight(ev, mc, daughters);
      } else if (ev.mc_type == 5) {
        a1_weight(ev, mc, daughters);
      }

      // check on negative weights
      if (ev.mc_weight_plus < 0.0 || ev.mc_weight_minus < 0.0) {
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
// helper to fill P4 and add energy & charge
void fill_collection (myEvent &ev, const auto &input_particles,
                            RVec<TLorentzVector> &out_p4) {

  out_p4.reserve(input_particles.size());
  
  for (const auto &p : input_particles) {
    ev.m_RecoCharge += p.charge;
    ev.m_RecoEnergy += p.energy;

    TLorentzVector tlv;
    tlv.SetPxPyPzE(p.momentum.x, p.momentum.y, p.momentum.z, p.energy);
    out_p4.push_back(tlv);
  }
}
// ==========================================
// helper
RVec<int> get_idx(const bool &hemi, const RVec<double> &costheta, const RVec<int> &ids){
  
  RVec<int> out;
  if(hemi){
    for(int i=0;i<costheta.size();i++){
      if (costheta[i]>=0) out.push_back(ids[i]);
    }
  } else if (!hemi) {
    for(int i=0;i<costheta.size();i++){
      if (costheta[i]<0) out.push_back(ids[i]);
    }
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
    //else if (p == 211 || p == 321 || p == 323)
    else if (p == 211)
      n_pi++; 
    else if (p == 111)
      n_pi0++;
    else if (p == 22)
      n_ph++;
  }
  // classification
  if (n_mu == 1 )
    return 1; // mu
  else if (n_el == 1 )
    return 2; // el
  else if (n_pi == 1 && n_pi0 == 0 && n_ph == 0)
    return 3; // pi
  else if (n_pi == 1 && n_pi0 == 1 && n_ph == 0)
    return 4; // rho
  //else if (n_pi == 1 && n_pi0 == 0 && n_ph == 2)
    //return 4; // rho (not a pi0 in the decay but directly gammas)
  else if (n_pi == 1 && n_pi0 == 2 && n_ph == 0)
    return 5; // a1 (1prong)
  else if (n_pi == 3)
    return 5; // a1 (3prong)

  return 0;
}

// ==========================================
int classify_lep(myEvent &ev, const RVec<int> &mu_idx, const RVec<int> &el_idx, const RVec<int> &rp2mc_idx) {
  // Check MUON: 1 mu, 0 others
  if (ev.n_mu == 1 && ev.n_el == 0) {

    // check RP2MC and get corresponding 1prong true MC energy
    int id = mu_idx[0];
    ev.mc_RP2MC_id = rp2mc_idx[id];
    return 1; // Type 1: Muon
  }
  // Check ELECTRON: 1 el, 0 others
  if (ev.n_mu == 0 && ev.n_el == 1 ) {

    // check RP2MC and get corresponding 1prong true MC energy
    int id = el_idx[0];
    ev.mc_RP2MC_id = rp2mc_idx[id];
    return 2; // Type 2: Electron
  }

  return 0;
}

// ==========================================
int classify_pion(myEvent &ev, const RVec<int> &pi_idx, const RVec<int> &rp2mc_idx) {

  // invariant mass check is done in main function
  // classification
  if (ev.n_pi == 1) {
    if (ev.n_ph == 0) {
      // check RP2MC and get corresponding 1prong true MC energy
      int id = pi_idx[0];
      ev.mc_RP2MC_id = rp2mc_idx[id];
      return 3; // Type 3: Single Pion
  } else if (ev.n_ph >= 1 && ev.n_ph <= 2) {

      // flag: 0,2 GeV < mass < 1.4 GeV
      if(ev.m_RecoMass>1.4 || ev.m_RecoMass<0.2) {
        ev.m_debug_mass = 11;
      }
      return 4; // Type 4: Rho (pi + 1-2 gamma)
  } else if (ev.n_ph >= 3) {

      // flag: 0,6 < mass < 1.8 GeV
      if(ev.m_RecoMass>1.8 || ev.m_RecoMass<0.6) {
        ev.m_debug_mass = 11;
      }
      ev.m_debug = 2;
      return 5; // Type 5: a1 -> pi + 2pi0 -> pi + 4gamma
    }
  }
  if (ev.n_pi == 3) {
  	ev.m_debug = 3;

    // flag: 0,6 < mass < 1.8 GeV
      if(ev.m_RecoMass>1.8 || ev.m_RecoMass<0.6) {
        ev.m_debug_mass = 11;
      }
    return 5; // Type 5: a1 (3-prong mode)
  }
  
  return 0;
}
  
  

RVec<int> get_type_safe(const RVec<myEvent> &evs, const bool masscheck) {
  RVec<int> out;
  out.reserve(evs.size());
  for (const auto &e : evs) {
    if(masscheck && e.m_debug_mass != 1 && e.m_debug_mass != 11){
      
      // impose 2 gev energy cut for pion decays
      if(e.m_type == 3 && e.m_RecoEnergy < 2.0){
        out.push_back(0);
        continue;
      }
      out.push_back(e.m_type);

    }
    else if (!masscheck && e.m_debug_mass != 1) out.push_back(e.m_type);
    else out.push_back(0);
  }
  return out;
}

// ==========================================
// RE-WEIGHTING FUNCTIONS
// ==========================================

double calc_Ptau(const TLorentzVector &p4_tau) {
  double costheta = p4_tau.CosTheta();
  double Ptau = -(SM_Atau * (1 + costheta * costheta) + 2 * SM_Atau * costheta) /
               (1 + costheta * costheta + 2 * SM_Atau * SM_Atau * costheta);
  return Ptau;
}

// Function to calculate cos(theta*) using Truth MC 4-vectors
double GetCosThetaStar(const TLorentzVector &p4_tau_lab,
                      const TLorentzVector &p4_pi_lab) {

  // BoostVector() returns beta of the particle (Lab -> Particle)
  TVector3 boost_to_rest = -p4_tau_lab.BoostVector();
  // copying and boosting
  TLorentzVector p4_pi_rest = p4_pi_lab;
  p4_pi_rest.Boost(boost_to_rest);

  // get tau's direction
  TVector3 tau_dir_lab = p4_tau_lab.Vect();

  // calculate angle between tau and pi boosted in tau rest frame
  double angle = p4_pi_rest.Vect().Angle(tau_dir_lab);

  return cos(angle);
}

void lepton_weight(myEvent &ev, const RVec<edm4hep::MCParticleData> &mc,
                 const RVec<int> &daughters) {

  double Ptau = 0.; // recalculated from tau p4
  const int tau_idx = ev.mc_tau_index;

  // use tau index to find tau directly
  if (tau_idx < 0 || tau_idx >= mc.size()) {
    cerr << "[ERROR]: Invalid tau index" << endl;
    return;
  }

  // tau found
  const auto &p = mc[tau_idx];
  TLorentzVector p4_tau_lab;
  p4_tau_lab.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
  Ptau = calc_Ptau(p4_tau_lab);
  ev.mc_Ptau = Ptau;
  // cycle daughters and store lepton true p4
  int pb = p.daughters_begin;
  int pe = p.daughters_end;
  TLorentzVector p4_lep_lab;
  for (int i = pb; i < pe; i++) {
    int dau_idx = daughters[i];
    const auto &dau = mc[dau_idx];
    if (abs(dau.PDG) == 13 || abs(dau.PDG) == 11 ) {
      // found lepton daughter
      ev.m_found = true;
      p4_lep_lab.SetXYZM(dau.momentum.x, dau.momentum.y, dau.momentum.z,
                        dau.mass);
    }
  }

  ev.mc_daughterP4 = p4_lep_lab;
  ev.mc_daughterMass = p4_lep_lab.M();

  double lep_energy = p4_lep_lab.E();
  // add photons if event is electron decay
  if(ev.mc_type == 2){
    for (int i = pb; i < pe; i++) {
      int dau_idx = daughters[i];
      const auto &dau = mc[dau_idx];
      if (abs(dau.PDG) == 22 ) {
        // found photon daughter
        TLorentzVector p4_photon;
        p4_photon.SetXYZM(dau.momentum.x, dau.momentum.y, dau.momentum.z,
                          dau.mass);
        lep_energy += p4_photon.E();
      }
    }
  }
  

  double x = lep_energy/E_TAU;

  double a = (5.0-9.0*x*x+4.0*x*x*x);
  double b = (1.0-9.0*x*x+8.0*x*x*x);
  // weight
  double w_plus = (1 + b/a) / (1 + Ptau * (b/a));
  double w_minus = (1 - b/a) / (1 + Ptau * (b/a));

  ev.mc_weight_plus = w_plus;
  ev.mc_weight_minus = w_minus;
  ev.mc_omega = x;
}

void pion_weight(myEvent &ev, const RVec<edm4hep::MCParticleData> &mc,
                 const RVec<int> &daughters) {

  const double charge = ev.m_RecoCharge;
  const int tau_idx = ev.mc_tau_index;
  double z = 0.; // costheta star
  double alpha = 1.;
  double Ptau = 0.; // recalculated from tau p4

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
  ev.mc_Ptau = Ptau;
  TLorentzVector p4_pi_lab;
  
  for (int i = pb; i < pe; i++) {
    int dau_idx = daughters[i];
    const auto &dau = mc[dau_idx];
    
    if (abs(dau.PDG) == 211) {
      // found pion daughter
      ev.m_found = true;
      
      TLorentzVector p_temp;
      p_temp.SetXYZM(dau.momentum.x, dau.momentum.y, dau.momentum.z,
                        dau.mass);
      p4_pi_lab += p_temp;
    }
    
  }
  ev.mc_daughterP4 = p4_pi_lab;
  ev.mc_daughterMass = p4_pi_lab.M();
  z = GetCosThetaStar(p4_tau_lab, p4_pi_lab);
  // weight
  double w_plus = (1 + alpha * z) / (1 + alpha * Ptau * z);
  double w_minus = (1 - alpha * z) / (1 + alpha * Ptau * z);

  ev.mc_weight_plus = w_plus;
  ev.mc_weight_minus = w_minus;
  ev.mc_omega = p4_pi_lab.E()/E_TAU;
}

void rho_weight(myEvent &ev, const RVec<edm4hep::MCParticleData> &mc,
                const RVec<int> &daughters) {

  const double charge = ev.m_RecoCharge;
  const int tau_idx = ev.mc_tau_index;
  double z = 0.;       // costheta star
  double alpha = 0.46; // recalculate
  double Ptau = 0.;

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
  ev.mc_Ptau = Ptau;
  TLorentzVector p4_rho_lab, p4_pi_lab, p4_pi0_lab;

  for (int i = pb; i < pe; i++) {
    // loop on daughters, find resonance, store p4
    int dau_idx = daughters[i];
    const auto &dau = mc[dau_idx];
    // save pion
    if (abs(dau.PDG) == 211) {
      ev.m_found = true;
      p4_pi_lab.SetXYZM(dau.momentum.x, dau.momentum.y, dau.momentum.z,
                        dau.mass);
      ev.mc_piP4 = p4_pi_lab;
    } else if (abs(dau.PDG) == 111) {
      // save pi0
      p4_pi0_lab.SetXYZM(dau.momentum.x, dau.momentum.y, dau.momentum.z,
                         dau.mass);
      ev.mc_pi0P4 = p4_pi0_lab;
    } 
  }
  p4_rho_lab = p4_pi_lab + p4_pi0_lab;
  ev.mc_daughterP4 = p4_rho_lab;
  double mRho = p4_rho_lab.M();
  ev.mc_daughterMass = mRho;
  z = GetCosThetaStar(p4_tau_lab, p4_rho_lab);
  alpha = (SM_TAU * SM_TAU - 2 * mRho * mRho) / (SM_TAU * SM_TAU + 2 * mRho * mRho);
  // weights
  ev.mc_weight_plus = (1 + alpha * z) / (1 + alpha * Ptau * z);
  ev.mc_weight_minus = (1 - alpha * z) / (1 + alpha * Ptau * z);
  ev.mc_omega = geometric_omega_rho(ev,p4_tau_lab,p4_rho_lab,p4_pi_lab,p4_pi0_lab);
}

void a1_weight(myEvent &ev, const RVec<edm4hep::MCParticleData> &mc,
               const RVec<int> &daughters) {

  const double charge = ev.m_RecoCharge;
  const int tau_idx = ev.mc_tau_index;
  double z = 0.;       // costheta star
  double alpha = 0.12; // recalculate ?? worth it??
  double Ptau = 0.;

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
  ev.mc_Ptau = Ptau;
  TLorentzVector p4_a1_lab;

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
      p4_a1_lab += p4_pi_lab;
    } else if (abs(dau.PDG) == 111) {
      // save pi0
      TLorentzVector p4_pi0_lab;
      p4_pi0_lab.SetXYZM(dau.momentum.x, dau.momentum.y, dau.momentum.z,
                         dau.mass);
      p4_a1_lab += p4_pi0_lab;
    }
  }
  ev.mc_daughterP4 = p4_a1_lab;
  ev.mc_daughterMass = p4_a1_lab.M();
  z = GetCosThetaStar(p4_tau_lab, p4_a1_lab);
  // weights
  ev.mc_weight_plus = (1 + alpha * z) / (1 + alpha * Ptau * z);
  ev.mc_weight_minus = (1 - alpha * z) / (1 + alpha * Ptau * z);
}

// optimal variable for rho decays
double calculate_omega_rho(myEvent &ev, TLorentzVector &p4_tau, 
                     const TLorentzVector &p4_rho, 
                     const TLorentzVector &p4_pip, 
                     const TLorentzVector &p4_pi0) {

  
  // get mass and energies from p4
  double m_rho = p4_rho.M();
  double E_rho = p4_rho.E();
  double E_tau = E_TAU;
  // double E_tau = p4_tau.E();
  double P_rho = p4_rho.P();

  // SAFETY CHECK
  if (m_rho < 2.0 * SM_PI + 0.001) {
     return -999.0; 
  }

  // angle of rho in tau rest frame
  // cos_psi_tau = (2x - 1 - m_rho^2/m_tau^2) / (1 - m_rho^2/m_tau^2)
  // where x = E_rho / E_tau
  
  double cos_psi_tau = (2.0 * E_rho/E_tau - 1.0 - (m_rho*m_rho)/(SM_TAU*SM_TAU) ) / (1.0 - (m_rho*m_rho)/(SM_TAU*SM_TAU) );

  if (cos_psi_tau > 1.0) {
    cos_psi_tau = 1.0;
    //flag bad event
    ev.m_debug = 99;
  }
  if (cos_psi_tau < -1.0) {
    cos_psi_tau = -1.0;
    //flag bad event
    ev.m_debug = 99;
  }

  // angle of charged pion in rho rest frame
  // cos_psi_rho = (m_rho / sqrt(m_rho^2 - 4m_pi^2)) * (E_pi_charged - E_pi_neutral) / P_rho
  // check sqrt argument
  if(m_rho*m_rho - 4.0*SM_PI*SM_PI <= 0){
    ev.m_debug = 99;
    // should check if these are events with missing photons or just wrong ids
    return -999;
  }
  double cos_psi_rho = (m_rho / sqrt(m_rho*m_rho - 4.0*SM_PI*SM_PI)) * (p4_pip.E() - p4_pi0.E()) / P_rho;

  if (cos_psi_rho > 1.0) {
    cos_psi_rho = 1.0;
    //flag bad event
    ev.m_debug = 99;
  }
  if (cos_psi_rho < -1.0) {
    cos_psi_rho = -1.0;
    //flag bad event
    ev.m_debug = 99;
  }

  if(cos_psi_rho > 1.0 || cos_psi_rho < -1.0){
    cerr << "[WARNING]: cos(psi_rho) out of bounds: " << cos_psi_rho << " , reco evt = " << ev.m_type << endl;
  }
  if(cos_psi_tau > 1.0 || cos_psi_tau < -1.0){
    cerr << "[WARNING]: cos(psi_tau) out of bounds: " << cos_psi_tau << " , reco evt = " << ev.m_type << endl;
  } 

  // compute Omega
  double psi_tau = acos(cos_psi_tau);

  // Wigner Rotation Angle eta 
  // tan(eta/2) = (m_rho/m_tau) * tan(psi_tau/2)
  double tan_theta_2 = tan(psi_tau / 2.0);
  double eta = 2.0 * atan( (m_rho / SM_TAU) * tan_theta_2 );

  // Precompute trig terms for w functions
  double cos_eta     = cos(eta);
  double sin_eta     = sin(eta);
  double cos_theta_2 = cos(psi_tau / 2.0); 
  double sin_theta_2 = sin(psi_tau / 2.0);

  // Calculate w_0 and w_1 components 
  // w0+
  double term0p = SM_TAU * cos_eta * cos_theta_2 + m_rho * sin_eta * sin_theta_2;
  double w0_plus = term0p * term0p;

  // w0-
  double term0m = SM_TAU * cos_eta * sin_theta_2 - m_rho * sin_eta * cos_theta_2;
  double w0_minus = term0m * term0m;

  // w1 +
  double term1p = SM_TAU * sin_eta * cos_theta_2 - m_rho* cos_eta * sin_theta_2;
  double w1_plus = (term1p * term1p) + (m_rho * m_rho * sin_theta_2 * sin_theta_2);

  // w1 -
  double term1m = SM_TAU * sin_eta * sin_theta_2 + m_rho * cos_eta * cos_theta_2;
  double w1_minus = (term1m * term1m) + (m_rho * m_rho * cos_theta_2 * cos_theta_2);

  // Calculate h functions 
  double h0 = 2.0 * cos_psi_rho * cos_psi_rho;
  double h1 = 1.0 - cos_psi_rho * cos_psi_rho;

  // Calculate W+ and W- 
  double W_plus  = w1_plus * h1 + w0_plus * h0 ;
  double W_minus = w1_minus * h1 + w0_minus * h0 ;


  double omega =  (W_plus - W_minus) / (W_plus + W_minus);

  if(omega > 1.0) cerr << "[WARNING]: omega > 1.0 (" << omega << ")" << endl;
  if(omega < -1.0) cerr << "[WARNING]: omega < -1.0 (" << omega << ")" << endl;

  return omega;
}

// omega_rho with p4 angles instead of kinematic variables, possible only for MC
// how to implement p4_tau for reco???? TODO
double geometric_omega_rho(myEvent &ev, TLorentzVector &p4_tau, 
                     const TLorentzVector &p4_rho, 
                     const TLorentzVector &p4_pip, 
                     const TLorentzVector &p4_pi0) {

  
  // get MC variables
  double m_rho = p4_rho.M();
  double E_rho = p4_rho.E();
  double E_tau = p4_tau.E();
  double P_rho = p4_rho.P();

  // angle of rho in tau rest frame
  
  double cos_psi_tau = GetCosThetaStar(p4_tau, p4_rho);

  if (cos_psi_tau > 1.0)  cos_psi_tau = 1.0;
  if (cos_psi_tau < -1.0) cos_psi_tau = -1.0;

  // aAngle of charged pion in rho rest frame
  
  double cos_psi_rho = GetCosThetaStar(p4_rho, p4_pip);

  if (cos_psi_rho > 1.0)  cos_psi_rho = 1.0;
  if (cos_psi_rho < -1.0) cos_psi_rho = -1.0;

  if(cos_psi_rho > 1.0 || cos_psi_rho < -1.0){
    cerr << "[WARNING]: cos(psi_rho) out of bounds: " << cos_psi_rho << " , reco evt = " << ev.m_type << endl;
  }
  if(cos_psi_tau > 1.0 || cos_psi_tau < -1.0){
    cerr << "[WARNING]: cos(psi_tau) out of bounds: " << cos_psi_tau << " , reco evt = " << ev.m_type << endl;
  } 

  double psi_tau = acos(cos_psi_tau);
  double tan_theta_2 = tan(psi_tau / 2.0);
  double eta = 2.0 * atan( (m_rho / SM_TAU) * tan_theta_2 );

  // Precompute trig terms for w functions
  double cos_eta     = cos(eta);
  double sin_eta     = sin(eta);
  double cos_theta_2 = cos(psi_tau / 2.0); 
  double sin_theta_2 = sin(psi_tau / 2.0);

  // Calculate w_0 and w_1 components 
  // w0+
  double term0p = SM_TAU * cos_eta * cos_theta_2 + m_rho * sin_eta * sin_theta_2;
  double w0_plus = term0p * term0p;

  // w0-
  double term0m = SM_TAU * cos_eta * sin_theta_2 - m_rho * sin_eta * cos_theta_2;
  double w0_minus = term0m * term0m;

  // w1 +
  double term1p = SM_TAU * sin_eta * cos_theta_2 - m_rho* cos_eta * sin_theta_2;
  double w1_plus = (term1p * term1p) + (m_rho * m_rho * sin_theta_2 * sin_theta_2);

  // w1 -
  double term1m = SM_TAU * sin_eta * sin_theta_2 + m_rho * cos_eta * cos_theta_2;
  double w1_minus = (term1m * term1m) + (m_rho * m_rho * cos_theta_2 * cos_theta_2);

  // Calculate h functions 
  double h0 = 2.0 * cos_psi_rho * cos_psi_rho;
  double h1 = 1.0 - cos_psi_rho * cos_psi_rho;

  // Calculate W+ and W- 
  double W_plus  = w1_plus * h1 + w0_plus * h0 ;
  double W_minus = w1_minus * h1 + w0_minus * h0 ;


  double omega =  (W_plus - W_minus) / (W_plus + W_minus);

  if(omega > 1.0) cerr << "[WARNING]: omega > 1.0 (" << omega << ")" << endl;
  if(omega < -1.0) cerr << "[WARNING]: omega < -1.0 (" << omega << ")" << endl;

  return omega;
}

void new_rho_weight(myEvent &ev, const RVec<edm4hep::MCParticleData> &mc,
                const RVec<int> &daughters) {

  const int tau_idx = ev.mc_tau_index;
  double Ptau = 0.;
  
  // tlv to be filled
  TLorentzVector p4_tau_lab, p4_rho_lab, p4_pip_lab, p4_pi0_lab;
  bool found_pip = false;
  bool found_pi0 = false;

  // Find Tau and Ptau
  if (tau_idx < 0 || tau_idx >= mc.size()) {
    cerr << "[ERROR]: Invalid tau index" << endl;
    return;
  }
  const auto &p = mc[tau_idx];
  p4_tau_lab.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
  
  Ptau = calc_Ptau(p4_tau_lab);
  ev.mc_Ptau = Ptau;

  // Cycle daughters to build Rho, Pi+, Pi0
  int pb = p.daughters_begin;
  int pe = p.daughters_end;

  for (int i = pb; i < pe; i++) {
    int dau_idx = daughters[i];
    const auto &dau = mc[dau_idx];
    
    // Create temp 4-vector
    TLorentzVector p4_dau;
    p4_dau.SetXYZM(dau.momentum.x, dau.momentum.y, dau.momentum.z, dau.mass);

    if (abs(dau.PDG) == 211 ) { 
      p4_pip_lab = p4_dau;
      ev.mc_piP4 = p4_dau;
      found_pip = true;
    } 
    // Neutral Pion (111)
    else if (abs(dau.PDG) == 111) {
      p4_pi0_lab = p4_dau;
      ev.mc_pi0P4 = p4_dau;
      found_pi0 = true;
    }
  }

  
  if (found_pip && found_pi0) {
    ev.m_found = true;
    
    // Reconstruct Rho 4-vector from daughters
    p4_rho_lab = p4_pip_lab + p4_pi0_lab;

    // Store MC info in event struct
    ev.mc_daughterP4 = p4_rho_lab;
    ev.mc_daughterMass = p4_rho_lab.M();

    double omega = calculate_omega_rho(ev,p4_tau_lab, p4_rho_lab, p4_pip_lab, p4_pi0_lab);

    // ---------------------------------------------------------
    // ASSIGN WEIGHTS
    // W ~ 1 + Ptau * omega
    // H=+1 (Ptau=+1) -> W ~ 1 + omega
    // H=-1 (Ptau=-1) -> W ~ 1 - omega
    // ---------------------------------------------------------
    
    
    double weight_denom = 1.0 + Ptau * omega;
    if (fabs(weight_denom) < 1e-6) weight_denom = 1.0; // safety

    ev.mc_weight_plus  = (1.0 + omega) / weight_denom;
    ev.mc_weight_minus = (1.0 - omega) / weight_denom;

  } else {
    // Fallback if daughters not found correctly
    ev.mc_weight_plus = 1.0;
    ev.mc_weight_minus = 1.0;
  }
}

// ==========================================
// EXTRACT VARIABLES
// ==========================================
RVec<double> get_lepton_e(const RVec<myEvent> &evs, const int mc_type,
                         const bool bool_mc, const int reco_type,
                         const bool bool_reco, const bool masscheck,
                          const bool asymmetric) {

  RVec<double> out;
  // choose only events with 1 hadronic tau + 1 leptonic tau if asked
  if(asymmetric){
    // skip any 'other' non-classified event
    if(evs[0].mc_type == 0 || evs[1].mc_type ==0) return out;
    // check 1 hadronic + 1 leptonic
    if( ( (evs[0].mc_type <=2) && (evs[1].mc_type <=2) ) ||
        ( (evs[0].mc_type >=3) && (evs[1].mc_type >=3) ) ) return out;
  }
  for (const auto &e : evs) {
    // check mc event
    if (bool_mc && e.mc_type != mc_type)
      continue;
    // check reco event
    if (bool_reco && e.m_type != reco_type)
      continue;
    // impose invariant mass check
    if (masscheck && e.m_debug_mass == 1)
      continue;
    // lepton energy
    if (e.n_mu > 0)
      out.push_back(e.m_muP4[0].E());
    else if (e.n_el > 0)
      out.push_back(e.m_elP4[0].E());
  }
  return out;
};
// ==========================================
RVec<double> get_hadron_e(const RVec<myEvent> &evs, const int mc_type, const bool bool_mc,
                          const int reco_type, const bool masscheck) {

  RVec<double> out;

  for (const auto &e : evs) {
  
    // impose invariant mass check
    if ( masscheck && (e.m_debug_mass == 1 || e.m_debug_mass == 11) ) 
      continue;

    // if bool_mc is false skip mc type check
    if (!bool_mc){
      if (e.m_type == reco_type) {
        out.push_back(e.m_RecoEnergy);
      }
      continue;
    }
    // from now we assume bool_mc is true
    
    // if reco && mc type are >0 collect energy where reco == mc
    if ( (mc_type >= 0) && (reco_type >= 0) ) {
      if (e.mc_type == mc_type && e.m_type == reco_type) {
        out.push_back(e.m_RecoEnergy);
      }
    } 
    // if reco && mc type are <0 collect any energy where reco != mc given a said reco decay
    else if ( (mc_type <= 0) && (reco_type <= 0) ) {
      if ( (e.m_type == reco_type) && (e.mc_type != e.m_type) ) {
        out.push_back(e.m_RecoEnergy);
      }
    }
  }

  return out;
};

// =========================================
RVec<double> get_optimal(RVec<myEvent> &evs, const int mc_type,
                         const bool bool_mc, const int reco_type,
                         const bool bool_reco, const bool masscheck){

  RVec<double> out; 

  for (auto &e : evs) {
    // check mc event
    if (bool_mc && e.mc_type != mc_type)
      continue;
    // check reco event
    if (bool_reco && e.m_type != reco_type)
      continue;
    // impose invariant mass check
    if ((masscheck && e.m_debug_mass == 1) || (masscheck && e.m_debug_mass == 11))
      continue;
    // optimal variable
    out.push_back(e.mc_omega);
    
  }

  return out;
}


// =========================================
RVec<double> get_reco_x(RVec<myEvent> &evs, const int mc_type,
                         const bool bool_mc, const int reco_type,
                         const bool bool_reco, const bool masscheck){

  RVec<double> out; 

  for (auto &e : evs) {
    // check mc event
    if (bool_mc && e.mc_type != mc_type)
      continue;
    // check reco event
    if (bool_reco && e.m_type != reco_type)
      continue;
    // impose invariant mass check
    if ((masscheck && e.m_debug_mass == 1) || (masscheck && e.m_debug_mass == 11))
      continue;
    // optimal variable
    out.push_back(e.m_RecoEnergy/E_TAU);
    
  }

  return out;
}

// ==========================================
RVec<double> get_omega_rho(RVec<myEvent> &evs, const int mc_type,
                         const bool bool_mc, const int reco_type,
                         const bool bool_reco, const bool masscheck,
                          const bool asymmetric) {

  RVec<double> out;
  // choose only events with 1 hadronic tau + 1 leptonic tau if asked
  if(asymmetric){
    // skip any 'other' non-classified event
    if(evs[0].mc_type == 0 || evs[1].mc_type ==0) return out;
    // check 1 hadronic + 1 leptonic
    if( ( (evs[0].mc_type <=2) && (evs[1].mc_type <=2) ) ||
        ( (evs[0].mc_type >=3) && (evs[1].mc_type >=3) ) ) return out;
  }
  for (auto &e : evs) {
    // check mc event
    if (bool_mc && e.mc_type != mc_type)
      continue;
    // check reco event
    if (bool_reco && e.m_type != reco_type)
      continue;
    // impose invariant mass check
    if ((masscheck && e.m_debug_mass == 1) || (masscheck && e.m_debug_mass == 11))
      continue;
    // impose cos limits with m_debug
    if (masscheck && e.m_debug == 99)
      continue;
    // omega_rho
    TLorentzVector p4_tau, p4_rho, p4_pip, p4_pi0;

    if (bool_reco || (bool_reco && bool_mc)){
      // fill with reco p4
      for (const auto &p : e.m_piP4) {
        p4_pip += p;
      }
      for (const auto &p : e.m_phP4) {
        p4_pi0 += p;
      }

      p4_rho = p4_pip + p4_pi0;
      p4_tau = e.mc_tauP4;
      double omega = calculate_omega_rho(e,p4_tau, p4_rho, p4_pip, p4_pi0);
      out.push_back(omega);
    }

    if (bool_mc && !bool_reco){
      // fill with mc p4
      p4_pip = e.mc_piP4;
      p4_pi0 = e.mc_pi0P4;
      p4_tau = e.mc_tauP4;
      p4_rho = e.mc_daughterP4;
      double omega = geometric_omega_rho(e,p4_tau, p4_rho, p4_pip, p4_pi0);
      out.push_back(omega);
    }
    
    
  }
  return out;
};


// ==========================================
RVec<double> get_photon_e(const RVec<myEvent> &evs, const int mc_type,
                         const bool bool_mc, const int reco_type,
                         const bool bool_reco, const bool masscheck,
                          const bool asymmetric) {

  RVec<double> out;
  // choose only events with 1 hadronic tau + 1 leptonic tau if asked
  if(asymmetric){
    // skip any 'other' non-classified event
    if(evs[0].mc_type == 0 || evs[1].mc_type ==0) return out;
    // check 1 hadronic + 1 leptonic
    if( ( (evs[0].mc_type <=2) && (evs[1].mc_type <=2) ) ||
        ( (evs[0].mc_type >=3) && (evs[1].mc_type >=3) ) ) return out;
  }
  for (const auto &e : evs) {
    // check mc event
    if (bool_mc && e.mc_type != mc_type)
      continue;
    // check reco event
    if (bool_reco && e.m_type != reco_type)
      continue;
    // impose invarian mass check
    if (masscheck && e.m_debug_mass == 1)
      continue;
    // photon energies
    if (e.n_ph > 0) {
      for (const auto &p4 : e.m_phP4) {
        out.push_back(p4.E());
      }
    }
  }
  return out;
};
// ==========================================
RVec<double> get_dressed_e(const RVec<myEvent> &evs, const int mc_type,
                         const bool bool_mc, const int reco_type,
                         const bool bool_reco, const bool masscheck,
                          const bool asymmetric) {

  RVec<double> out;
  // choose only events with 1 hadronic tau + 1 leptonic tau if asked
  if(asymmetric){
    // skip any 'other' non-classified event
    if(evs[0].mc_type == 0 || evs[1].mc_type ==0) return out;
    // check 1 hadronic + 1 leptonic
    if( ( (evs[0].mc_type <=2) && (evs[1].mc_type <=2) ) ||
        ( (evs[0].mc_type >=3) && (evs[1].mc_type >=3) ) ) return out;
  }
  for (const auto &e : evs) {
    // check mc event
    if (bool_mc && e.mc_type != mc_type)
      continue;
    // check reco event
    if (bool_reco && e.m_type != reco_type)
      continue;
    // impose invariant mass check
    if (masscheck && e.m_debug_mass == 1)
      continue;
    // collect dressed energy (if the event is mu or el the other P4 vectors are empty)
    double temp=0.0;
    for(const auto &p4 : e.m_elP4) {
      temp += p4.E();
    }
    for(const auto &p4 : e.m_muP4) {
      temp += p4.E();
    }
    for(const auto &p4 : e.m_phP4) {
      temp += p4.E();
    }
    out.push_back(temp);
  }

  return out;
};
// ==========================================
RVec<double> get_rp2mc_e(const RVec<myEvent> &evs, const int mc_type,
                         const bool bool_mc, const int reco_type,
                         const bool bool_reco, const bool masscheck,
                          const bool asymmetric) {

  RVec<double> out;
  // choose only events with 1 hadronic tau + 1 leptonic tau if asked
  if(asymmetric){
    // skip any 'other' non-classified event
    if(evs[0].mc_type == 0 || evs[1].mc_type ==0) return out;
    // check 1 hadronic + 1 leptonic
    if( ( (evs[0].mc_type <=2) && (evs[1].mc_type <=2) ) ||
        ( (evs[0].mc_type >=3) && (evs[1].mc_type >=3) ) ) return out;
  }
  for (const auto &e : evs) {
    // check mc event
    if (bool_mc && e.mc_type != mc_type)
      continue;
    // check reco event
    if (bool_reco && e.m_type != reco_type)
      continue;
    // impose invariant mass check
    if (masscheck && e.m_debug_mass == 1)
      continue;
    // return rp2mc_e
    out.push_back(e.mc_RP2MC_e);
  }

  return out;
};
// ==========================================
RVec<double> get_weights(const int sign, const RVec<myEvent> &evs,
                        const int mc_type, const bool bool_mc,
                        const int reco_type, const bool bool_reco, const bool masscheck,
                        const bool asymmetric) {

  RVec<double> out;
  // choose only events with 1 hadronic tau + 1 leptonic tau if asked
  if(asymmetric){
    // skip any 'other' non-classified event
    if(evs[0].mc_type == 0 || evs[1].mc_type ==0) return out;
    // check 1 hadronic + 1 leptonic
    if( ( (evs[0].mc_type <=2) && (evs[1].mc_type <=2) ) ||
        ( (evs[0].mc_type >=3) && (evs[1].mc_type >=3) ) ) return out;
  }
  for (const auto &e : evs) {
    // check mc event
    if (bool_mc && e.mc_type != mc_type)
      continue;
    // check reco event
    if (bool_reco && e.m_type != reco_type)
      continue;
    // impose invariant mass check
    if ((masscheck && e.m_debug_mass == 1) || (masscheck && e.m_debug_mass == 11))
      continue;
    // get weight based on sign passed
    if (sign > 0)
      out.push_back(e.mc_weight_plus);
    else if (sign < 0)
      out.push_back(e.mc_weight_minus);
  }
  return out;
};
// ==========================================
RVec<double> get_invariant_mass(const RVec<myEvent> &evs, const int mc_type,
                               const bool bool_mc, const int reco_type,
                               const bool bool_reco) {

  RVec<double> out;
  for (const auto &e : evs) {
    // check mc event
    if (bool_mc && e.mc_type != mc_type)
      continue;
    // check reco event
    if (bool_reco && e.m_type != reco_type)
      continue;
    // invariant mass
    // check debug
    if(e.m_debug_mass == 1 || e.m_debug_mass == 11){
      continue;
    }
    out.push_back(e.m_RecoMass);
  }
  return out;
};
// ==========================================
RVec<double> get_mass_pull(const RVec<myEvent> &evs, const int mc_type,
                       const bool bool_mc, const int reco_type,
                       const bool bool_reco) {

  RVec<double> out;
  for (const auto &e : evs) {
    // check mc event
    if (bool_mc && e.mc_type != mc_type)
      continue;
    // check reco event
    if (bool_reco && e.m_type != reco_type)
      continue;
    // check mass debug
    if(e.m_debug_mass == 1 || e.m_debug_mass == 11){
      continue;
    }
    // meson mass
    out.push_back(e.mc_daughterMass - e.m_RecoMass);
  }
  return out;
};
// ==========================================
RVec<double> get_MCdaughter_e(const RVec<myEvent> &evs, const int mc_type,
                       const bool bool_mc, const int reco_type,
                       const bool bool_reco, const bool masscheck,
                        const bool asymmetric) {

  RVec<double> out;
  // choose only events with 1 hadronic tau + 1 leptonic tau if asked
  if(asymmetric){
    // skip any 'other' non-classified event
    if(evs[0].mc_type == 0 || evs[1].mc_type ==0) return out;
    // check 1 hadronic + 1 leptonic
    if( ( (evs[0].mc_type <=2) && (evs[1].mc_type <=2) ) ||
        ( (evs[0].mc_type >=3) && (evs[1].mc_type >=3) ) ) return out;
  }
  for (const auto &e : evs) {
    // check mc event
    if (bool_mc && e.mc_type != mc_type)
      continue;
    // check reco event
    if (bool_reco && e.m_type != reco_type)
      continue;
    // impose invariant mass check
    if (masscheck && e.m_debug_mass == 1)
      continue;
    // true daughter MC energies
    out.push_back(e.mc_daughterP4.E());
  }
  return out;
};
// ==========================================
RVec<double> get_MCdaughter_x(const RVec<myEvent> &evs, const int mc_type,
                       const bool bool_mc, const int reco_type,
                       const bool bool_reco, const bool masscheck,
                        const bool asymmetric) {

  RVec<double> out;
  // choose only events with 1 hadronic tau + 1 leptonic tau if asked
  if(asymmetric){
    // skip any 'other' non-classified event
    if(evs[0].mc_type == 0 || evs[1].mc_type ==0) return out;
    // check 1 hadronic + 1 leptonic
    if( ( (evs[0].mc_type <=2) && (evs[1].mc_type <=2) ) ||
        ( (evs[0].mc_type >=3) && (evs[1].mc_type >=3) ) ) return out;
  }
  for (const auto &e : evs) {
    // check mc event
    if (bool_mc && e.mc_type != mc_type)
      continue;
    // check reco event
    if (bool_reco && e.m_type != reco_type)
      continue;
    // impose invariant mass check
    if (masscheck && e.m_debug_mass == 1)
      continue;
    // true daughter MC energies
    out.push_back(e.mc_daughterP4.E()/e.mc_tauP4.E());
  }
  return out;
};
// ==========================================
RVec<double> get_MCdaughter_mass(const RVec<myEvent> &evs, const int mc_type,
                       const bool bool_mc, const int reco_type,
                       const bool bool_reco) {

  RVec<double> out;
  for (const auto &e : evs) {
    // check mc event
    if (bool_mc && e.mc_type != mc_type)
      continue;
    // check reco event
    if (bool_reco && e.m_type != reco_type)
      continue;
    // meson mass
    out.push_back(e.mc_daughterMass);
  }
  return out;
};
// ==========================================
// DEBUG
// ==========================================

// ==========================================
RVec<double> get_Ptau(const RVec<myEvent> &evs,
                        const int mc_type, const bool bool_mc,
                        const int reco_type, const bool bool_reco, const bool masscheck,
                        const bool asymmetric) {

  RVec<double> out;
  // choose only events with 1 hadronic tau + 1 leptonic tau if asked
  if(asymmetric){
    // skip any 'other' non-classified event
    if(evs[0].mc_type == 0 || evs[1].mc_type ==0) return out;
    // check 1 hadronic + 1 leptonic
    if( ( (evs[0].mc_type <=2) && (evs[1].mc_type <=2) ) ||
        ( (evs[0].mc_type >=3) && (evs[1].mc_type >=3) ) ) return out;
  }
  for (const auto &e : evs) {
    // check mc event
    if (bool_mc && e.mc_type != mc_type)
      continue;
    // check reco event
    if (bool_reco && e.m_type != reco_type)
      continue;
    // impose invariant mass check
    if (masscheck && e.m_debug_mass == 1)
      continue;
    
    out.push_back(e.mc_Ptau);
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
        if (e.m_piP4[0].E() > 2.0)
          mask.push_back(1);
        else
          mask.push_back(0);
      }
    }
    // pass all other events
    else {
      mask.push_back(1);
    }
  }
  return mask;
}

RVec<int> get_weight_mask(const RVec<myEvent> &evs) {
  RVec<int> mask;
  for (const auto &e : evs) {
    // apply filter to pion events
    if (e.m_type == 3) {
      // safety check
      if (e.m_piP4.empty()) {
        mask.push_back(0);
      } else {
        // weight check
        if (e.mc_weight_plus < 0 || e.mc_weight_minus < 0)
          mask.push_back(0);
        else
          mask.push_back(1);
      }
    }
    // pass all other events
    else {
      mask.push_back(1);
    }
  }
  return mask;
}

RVec<int> get_debug_mass(const RVec<myEvent> &evs, const int reco_type) {
  RVec<int> out;
  for (const auto &e : evs) {
    // check reco event
    if (e.m_type != reco_type)
      continue;
    //check debug_mass
    if(e.m_debug_mass != 11)
      continue;
    out.push_back(e.mc_type);
  }
  return out;
}

RVec<int> get_type_debugmass(const RVec<myEvent> &evs,
                               const int mc_type, const bool bool_mc,
                               const int reco_type, const bool bool_reco) {
  RVec<int> out;
  for (const auto &e : evs) {
    // check mc event
    if (bool_mc && e.mc_type != mc_type)
      continue;
    // check reco event
    if (bool_reco && e.m_type != reco_type)
      continue;
    // check for invariant mass 
    if (e.m_debug_mass == 0)
      continue;
    // push back event type
    // note this can only work if we remove the reset of type in case of high mass
    // otherwise all types become 0
    out.push_back(e.m_type);
  }
  return out;
}

RVec<RVec<int>> get_debug_daughters(const RVec<myEvent> &evs,
                               const int mc_type, const bool bool_mc,
                               const int reco_type, const bool bool_reco) {
  RVec<RVec<int>> out;
  for (const auto &e : evs) {
    // check mc event
    if (bool_mc && e.mc_type != mc_type)
      continue;
    // check reco event (inverse)
    if (bool_reco && e.m_type == reco_type) // inverse search, look for events which are type A on montecarlo and NOT A in reco
      continue;
    // check for invariant mass 
    
    out.push_back(e.mc_daughters);
  }
  return out;
}

RVec<double> get_debug_costheta(const RVec<myEvent> &evs,
                               const int mc_type, const bool bool_mc,
                               const int reco_type, const bool bool_reco) {
  RVec<double> out;
  for (const auto &e : evs) {
    // check mc event
    if (bool_mc && e.mc_type != mc_type)
      continue;
    // check reco event 
    if (bool_reco && e.m_type != reco_type) 
      continue;

    out.push_back(e.thrust_costheta);
  }
  return out;
}

RVec<double> get_debug_phi(const RVec<myEvent> &evs,
                               const int mc_type, const bool bool_mc,
                               const int reco_type, const bool bool_reco) {
  RVec<double> out;
  for (const auto &e : evs) {
    // check mc event
    if (bool_mc && e.mc_type != mc_type)
      continue;
    // check reco event 
    if (bool_reco && e.m_type != reco_type) 
      continue;

    out.push_back(e.thrust_phi);
  }
  return out;
}


// ==========================================
RVec<debugInfo> get_debug_info(const RVec<myEvent> &evs, const int deb) {
  RVec<debugInfo> out;
  for (const auto &e : evs) {
    // filter by m_debug value
    if (e.m_debug == deb) {
      debugInfo info;
      info.m_debug = e.m_debug;
      info.m_type = e.m_type;
      info.mc_type = e.mc_type;
      out.push_back(info);
    }
  }
  return out;
}

// ==========================================
RVec<int> get_debug_info_mc_type(const RVec<debugInfo> &infos) {
  RVec<int> out;
  for (const auto &info : infos) {
    out.push_back(info.mc_type);
  }
  return out;
}

} // namespace Ztautau
