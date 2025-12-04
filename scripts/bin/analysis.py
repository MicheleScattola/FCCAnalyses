#Optional: List of analysis packages to load in runtime
analysesList = ['Ztautau']

#List of processes
processList = {
    'p8_ee_Ztautau_ecm91':{}
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
#prodTag     = "FCCee/winter2023/IDEA/"
# Link to the dictonary that contains all the cross section informations etc... (mandatory)
procDict = "FCCee_procDict_winter2023_IDEA.json"

#output directory
outputDir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/analysis/treemaker/"
outputName = "p8_ee_Ztautau_ecm91"
#input directory
inputDir    = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"

#Optional
nCPUS       = -1
#runBatch    = False
#batchQueue = "longlunch"
#compGroup = "group_u_FCC.local_gen"

#Optional test file , run with --test
testFile = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/p8_ee_Ztautau_ecm91/events_080694422.root"

#operations on the TTree
class RDFanalysis():
    #__________________________________________________________
    #analysers function to define the analysers to process
    def analysers(df):
        df2 = (
            df
                
        
				#MC truth and Reconstructed particles
				.Alias("Particle0", "Particle#0.index")
				.Alias("Particle1", "Particle#1.index")
				.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
				.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
				.Alias("rps", "ReconstructedParticles")
				.Define("rps_types","ReconstructedParticle::get_type(rps)")
				.Define("rps_MC_index",
  "ReconstructedParticle2MC::getRP2MC_index(MCRecoAssociations0, MCRecoAssociations1, rps)")
				
				#PDG ids of reconstructed particles
                .Define("AssociatedMCpdg","RVec<int> pdgs; for(auto idx : MCRecoAssociations1) pdgs.push_back(Particle[idx].PDG); return pdgs;")
                .Filter("rps.size()>2")
                
				
				#####
				# THRUST
				#####
				.Define("RP_e", "ReconstructedParticle::get_e(rps)")
				.Define("RP_px", "ReconstructedParticle::get_px(rps)")
				.Define("RP_py", "ReconstructedParticle::get_py(rps)")
				.Define("RP_pz", "ReconstructedParticle::get_pz(rps)")
				.Define("RP_charge", "ReconstructedParticle::get_charge(rps)")

				.Define("EVT_thrust0",      'Algorithms::minimize_thrust("Minuit2","Migrad")(RP_px, RP_py, RP_pz)')
				.Define("RP_thrustangle0",  'Algorithms::getAxisCosTheta(EVT_thrust0, RP_px, RP_py, RP_pz)')
				.Define("RP_thrustcostheta0","return EVT_thrust0[5]/ sqrt(EVT_thrust0[1]*EVT_thrust0[1] + EVT_thrust0[3]*EVT_thrust0[3] + EVT_thrust0[5]*EVT_thrust0[5])")
				.Define("RP_thrustphi0", "return atan2(EVT_thrust0[3],EVT_thrust0[1])")
				# calculating thrust with respect to charge in hemispheres
				.Define("EVT_thrust", "Ztautau::getThrustPointing(RP_charge, EVT_thrust0, RP_thrustangle0)")
				.Define("RP_thrustangle",   'Algorithms::getAxisCosTheta(EVT_thrust, RP_px, RP_py, RP_pz)')
				.Define("RP_thrustcostheta","return EVT_thrust[5]/ sqrt(EVT_thrust[1]*EVT_thrust[1]   + EVT_thrust[3]*EVT_thrust[3]   + EVT_thrust[5]*EVT_thrust[5])")
				.Define("RP_thrustphi", "return atan2(EVT_thrust[3],EVT_thrust[1])")
				
				#####
				# RECONSTRUCTED PARTICLES
				#####
				
				.Define("rps_pos", "ReconstructedParticle::sel_axis(true)(RP_thrustangle, rps)")
				.Define("rps_neg","ReconstructedParticle::sel_axis(false)(RP_thrustangle, rps)")
				.Define("Q_pos", 'Algorithms::getAxisCharge(true, 0.)(RP_thrustangle, RP_charge, RP_px, RP_py, RP_pz)')
				.Define("Q_neg", 'Algorithms::getAxisCharge(false,0.)(RP_thrustangle, RP_charge, RP_px, RP_py, RP_pz)')
				
				
                # Muons
				.Alias("Muon0", "Muon#0.index")
				.Define("muons",   "ReconstructedParticle::get(Muon0, rps)")
				.Define("muons_costheta","Ztautau::get_elements_by_index(RP_thrustangle,Muon0)")

				# Electrons
				.Alias("Electron0", "Electron#0.index")
				.Define("electrons",    "ReconstructedParticle::get(Electron0, rps)")
				.Define("electrons_costheta","Ztautau::get_elements_by_index(RP_thrustangle,Electron0)")

				# Photons
				.Alias("Photon0", "Photon#0.index")
				.Define("photons",   "ReconstructedParticle::get(Photon0, rps)")
				.Define("photons_e", "ReconstructedParticle::get_e(photons)")
				.Define("photons_costheta","Ztautau::get_elements_by_index(RP_thrustangle,Photon0)")
				
				# defining pions as charged hadrons with mass selection
				# selecting candidates (possibly mistaken with a K+ )
				.Define("pions_charged_ids", "Ztautau::sel_pions_id(rps,1)")
				.Define("pions_charged","ReconstructedParticle::get(pions_charged_ids, rps)")
				.Define("pi_charged_n","ReconstructedParticle::get_n(pions_charged)")
				.Define("pi_charged_e","ReconstructedParticle::get_e(pions_charged)")
				.Define("pi_costheta","Ztautau::get_elements_by_index(RP_thrustangle,pions_charged_ids)")
				
				#####
				# EVENTS IDENTIFICATION
				#####
				.Define("event_type_all","Ztautau::get_event( muons,muons_costheta, electrons, electrons_costheta, pions_charged, pi_costheta, photons, photons_costheta,10.)")
				.Define("event_type_reco","RVec<int> {event_type_all[0].type,event_type_all[1].type}")
				
				#####
				# MC IDENTIFICATION
				#####
				
				#TAU-
				.Define("TauMinus_el_ids",
				  "MCParticle::get_indices(15,  {11, -12, 16}, true, false, false, false)(Particle, Particle1)")

				.Define("TauMinus_mu_ids",
				  "MCParticle::get_indices(15,  {13, -14, 16}, true, false, false, false)(Particle, Particle1)")

				.Define("TauMinus_pi_ids",
				  "MCParticle::get_indices(15,  {-211, 16},     true, false, false, false)(Particle, Particle1)")

				.Define("TauMinus_rho_ids",
				  "MCParticle::get_indices(15,  {-211,111,16},  false, false, false, false)(Particle, Particle1)")
				
				#TAU+
				.Define("TauPlus_el_ids",
				  "MCParticle::get_indices(-15, {-11, 12, -16}, true, false, false, false)(Particle, Particle1)")

				.Define("TauPlus_mu_ids",
				  "MCParticle::get_indices(-15, {-13, 14, -16}, true, false, false, false)(Particle, Particle1)")

				.Define("TauPlus_pi_ids",
				  "MCParticle::get_indices(-15, {211, -16},     true, false, false, false)(Particle, Particle1)")

				.Define("TauPlus_rho_ids",
				  "MCParticle::get_indices(-15, {211,111,-16},  false, false, false, false)(Particle, Particle1)")	  
				
				.Define("MC_event_neg","Ztautau::MC_classified(TauMinus_mu_ids, TauMinus_el_ids, TauMinus_pi_ids, TauMinus_rho_ids)")
				.Define("MC_event_pos","Ztautau::MC_classified(TauPlus_mu_ids, TauPlus_el_ids, TauPlus_pi_ids, TauPlus_rho_ids)")
				.Define("MC_event","return RVec<int> {MC_event_pos, MC_event_neg}")
				
				#####
				# Pi0 resonance
				#####
			
				
				#####
				# Taus energy in collinear approximation
				#####
				
				.Define("missing_tlv","Ztautau::missingTLV(91.,rps,0.)")
				.Define("taus_e_collinear","Ztautau::collinear_approx(rps_pos,rps_neg,91.18)")
				.Define("tau_pos_e","taus_e_collinear[0]")
				.Define("tau_neg_e","taus_e_collinear[1]")
				
				#####
				# properties of MC photons from pi0 decay
				#####
				.Define("tau_to_rho_ids",
				  "MCParticle::get_indices(15,  {-211,111,16},  false, true, true, false)(Particle, Particle1)")
				.Define("pi0_id","tau_to_rho_ids[2]")
				.Define("pi0_to_gg_ids","MCParticle::get_indices_MotherByIndex(pi0_id,				 {22,22},true,true,false,Particle,Particle1)")
				.Define("gammagamma_rps",  "ReconstructedParticle2MC::selRP_matched_to_list( pi0_to_gg_ids, MCRecoAssociations0,MCRecoAssociations1,rps,Particle)")
				.Define("gammagamma_e","ReconstructedParticle::get_e(gammagamma_rps)")
				.Define("MC_gammagamma","Ztautau::getMC(pi0_to_gg_ids,Particle)")
				.Define("MC_gammagamma_e","MCParticle::get_e(MC_gammagamma)")
				
				# 1 pi + 1 gamma reconstructed, TrueMC vs ALL events
				.Define("pi1gamma1_true","Ztautau::study_ph(MC_event,event_type_all,true,1,false,4)")
				.Define("pi1gamma1_false","Ztautau::study_ph(MC_event,event_type_all,false,1,false,4)")
				# 1pi + 2 gamma reconstructed, TrueMC vs ALL events
				.Define("pi1gamma2_true","Ztautau::study_ph(MC_event,event_type_all,true,2,false,4)")
				.Define("pi1gamma2_false","Ztautau::study_ph(MC_event,event_type_all,false,2,false,4)")
				
				# NON-PI0 events which are reconstructed as 1pi + Ngammas
				.Define("wrongrho_1ph","Ztautau::study_ph(MC_event,event_type_all,false,1,true,4)")
				.Define("wrongrho_2ph","Ztautau::study_ph(MC_event,event_type_all,false,2,true,4)")
				# distributions for the sum of 2 gammas
				.Define("sum_wrong2ph","Ztautau::study_ph_sum(MC_event,event_type_all,false,2,true,4)")
				.Define("sum_true2ph","Ztautau::study_ph_sum(MC_event,event_type_all,false,2,false,4)")
				
				# search for parent pdgs in wrong non pi0 events in 1pi+1/2 gammas
				.Define("myPh","Ztautau::buildRPTruthCollection(rps,Particle,Particle0, rps_MC_index,RP_thrustangle,Photon0,MC_event,event_type_reco)")
				.Define("myPh_parPDG",
					"ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> out;"
					"out.reserve(myPh.size());"
					"for (const auto &p : myPh) out.push_back(p.mc_parent_pdg);"
					"return out;")
				  
                
				.Define("myPh_parGenStatus",
					"ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> out;"
					"out.reserve(myPh.size());"
					"for (const auto &p : myPh) out.push_back(p.mc_parent_genStatus);"
					"return out;")
				.Define("myPh_parSize",
					"ROOT::VecOps::RVec<int> out;"
					"out.reserve(myPh.size());"
					"for (const auto &p : myPh) out.push_back(p.ntot);"
					"return out;")
					
				.Define("testGen","MCParticle::get_genStatus(Particle)")
				.Define("testPDG","MCParticle::get_pdg(Particle)")
				
				# my particle collections
				.Define("myEl","Ztautau::buildRPTruthCollection(rps,Particle,Particle0, rps_MC_index,RP_thrustangle,Electron0,MC_event,event_type_reco)")
				.Define("myMu","Ztautau::buildRPTruthCollection(rps,Particle,Particle0, rps_MC_index,RP_thrustangle,Muon0,MC_event,event_type_reco)")
				.Define("myPi","Ztautau::buildRPTruthCollection(rps,Particle,Particle0, rps_MC_index,RP_thrustangle,pions_charged_ids,MC_event,event_type_reco)")
				
				.Define("parentPDG1","Ztautau::parentPDG(myPh,event_type_all,MC_event,1,4)")
				.Define("parentPDG2","Ztautau::parentPDG(myPh,event_type_all,MC_event,2,4)")
				
				# contributions to 1 gamma reconstructed ALL EVTS
				.Define("el1g_all","Ztautau::pdgtype_ph_e(myPh,event_type_all,MC_event,1,3,false,11,true)")
				.Define("mu1g_all","Ztautau::pdgtype_ph_e(myPh,event_type_all,MC_event,1,3,false,13,true)")
				.Define("tau1g_all","Ztautau::pdgtype_ph_e(myPh,event_type_all,MC_event,1,3,false,15,true)")
				.Define("ph1g_all","Ztautau::pdgtype_ph_e(myPh,event_type_all,MC_event,1,3,false,22,true)")
				.Define("pi01g_all","Ztautau::pdgtype_ph_e(myPh,event_type_all,MC_event,1,3,false,111,true)")
				
				# contributions to 1 gamma reconstructed TRUE
				.Define("el1g_true","Ztautau::pdgtype_ph_e(myPh,event_type_all,MC_event,1,3,true,11,true)")
				.Define("mu1g_true","Ztautau::pdgtype_ph_e(myPh,event_type_all,MC_event,1,3,true,13,true)")
				.Define("tau1g_true","Ztautau::pdgtype_ph_e(myPh,event_type_all,MC_event,1,3,true,15,true)")
				.Define("ph1g_true","Ztautau::pdgtype_ph_e(myPh,event_type_all,MC_event,1,3,true,22,true)")
				.Define("pi01g_true","Ztautau::pdgtype_ph_e(myPh,event_type_all,MC_event,1,3,true,111,true)")
				
				.Define("rhoaspi","Ztautau::confusion_ph_e(myPh,4,true,3,true)")
				.Define("otherasrho","Ztautau::confusion_ph_e(myPh,0,true,4,true)")
				.Define("rhoasother","Ztautau::confusion_ph_e(myPh,4,true,0,true)")
				
				
				
				.Filter("Q_pos>0 && Q_neg<0")
				  
                )
		
				        
        
        df2.Display(["myPh_parPDG","myPh_parGenStatus"],20).Print()
        df2.Display(["testPDG","testGen"],5).Print()
        df2.Display(["rps_types","RP_e","photons_e"],20).Print()

        #df2.Display(["MC_gammagamma_e","gammagamma_e","photons_e"],20).Print()
        #df2.Display(["tau_to_rho_ids","pi0_id","pi0_to_gg_ids"],20).Print()
        #print(df2.GetColumnType("event"))
        
        return df2
       

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
        	"RP_thrustangle0",
        	"RP_thrustangle",
        	"RP_thrustphi0",
        	"RP_thrustphi",
        	"RP_thrustcostheta0",
        	"RP_thrustcostheta",
        	"taus_e_collinear",
        	"tau_pos_e",
        	"tau_neg_e",
        	"gammagamma_e",
        	"photons_e",
        	"MC_gammagamma_e",
        	"MC_event",
        	"event_type_reco",
        	"pi1gamma1_true",
        	"pi1gamma1_false",
        	"pi1gamma2_true",
        	"pi1gamma2_false",
        	"wrongrho_1ph",
        	"wrongrho_2ph",
        	"sum_wrong2ph",
        	"sum_true2ph",
        	"parentPDG1",
        	"parentPDG2",
        	"el1g_all",
        	"mu1g_all",
        	"tau1g_all",
        	"ph1g_all",
        	"pi01g_all",
        	"el1g_true",
        	"mu1g_true",
        	"tau1g_true",
        	"ph1g_true",
        	"pi01g_true",
        	"rhoaspi",
        	"otherasrho",
        	"rhoasother"
        	
        	]
        return branchList

