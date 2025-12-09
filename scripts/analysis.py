#Optional: List of analysis packages to load in runtime
#analysesList = ['Ztautau']

#List of processes
processList = {
    'p8_ee_Ztautau_ecm91':{'fraction':0.02},
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
#prodTag     = "FCCee/winter2023/IDEA/"
# Link to the dictonary that contains all the cross section informations etc... (mandatory)
procDict = "FCCee_procDict_winter2023_IDEA.json"

#output directory
outputDir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/treemaker/"
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
                
                #initial basic cuts
                .Filter("rps.size()>=2")
				
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
				
				.Filter("Q_pos>0 && Q_neg<0")
				
				
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
				.Define("myEvent","Ztautau::myget_event(Muon0,Electron0,pions_charged_ids,Photon0,rps,RP_thrustangle,Particle,Particle1)")
                
				.Define("event_type_reco","Ztautau::get_type_safe(myEvent)")
				
				#####
				# MC IDENTIFICATION
				#####
                .Define("MC_event","RVec<int> {myEvent[0].m_MCtype,myEvent[1].m_MCtype}")
                .Define("debug","RVec<int> {myEvent[0].m_debug,myEvent[1].m_debug}")
                .Define("debug_mass","RVec<float> {myEvent[0].m_MCPtau,myEvent[1].m_MCPtau}")
                
				#.Define("weights_plus","RVec<float> {myEvent[0].m_MCweight_plus,myEvent[1].m_MCweight_plus}")
                #.Define("weights_minus","RVec<float> {myEvent[0].m_MCweight_minus,myEvent[1].m_MCweight_minus}")
                #.Define("found","RVec<int> {myEvent[0].m_found,myEvent[1].m_found}")
                
				# masks for later cuts and selections
				.Define("pi_mask", "Ztautau::get_pi_mask(myEvent)")
				.Define("weight_mask", "Ztautau::get_weight_mask(myEvent)")
				
				# parameters are: event collection, MC type, bool MC type, reco type, bool reco type
				# confusion matrix
				.Define("ElAsPi_e","Ztautau::get_hadron_e(myEvent,2,true,3,true)")
				.Define("MuAsPi_e","Ztautau::get_hadron_e(myEvent,1,true,3,true)")
				.Define("RhoAsPi_e","Ztautau::get_hadron_e(myEvent,4,true,3,true)")
				.Define("PiAsPi_e","Ztautau::get_hadron_e(myEvent,3,true,3,true)")
				.Define("A1AsPi_e","Ztautau::get_hadron_e(myEvent,5,true,3,true)")
				.Define("OtAsPi_e","Ztautau::get_hadron_e(myEvent,0,true,3,true)")
				# pi signal and weights
				#.Define("pi_sgn","Ztautau::get_MCpi_e(myEvent,3,true,3,false)/45.5")
				#.Define("w_plus","Ztautau::get_weights(1,myEvent,3,true,3,false)")
				#.Define("w_minus","Ztautau::get_weights(-1,myEvent,3,true,3,false)")
				.Define("pi_sgn","Ztautau::get_hadron_e(myEvent,3,false,3,true)/45.5")
				.Define("w_plus","Ztautau::get_weights(1,myEvent,3,false,3,true)")
				.Define("w_minus","Ztautau::get_weights(-1,myEvent,3,false,3,true)")
				# lepton signals
				.Define("el_sgn","Ztautau::get_lepton_e(myEvent,2,false,2,true)/45.5")
				.Define("mu_sgn","Ztautau::get_lepton_e(myEvent,1,false,1,true)/45.5")
				# a1 confusion
				#.Define("pi_m","Ztautau::get_invariant_mass(myEvent,3,false,3,true)")
				#.Define("rho_m","Ztautau::get_invariant_mass(myEvent,4,false,4,true)")
				#.Define("a1_m","Ztautau::get_invariant_mass(myEvent,5,false,5,true)")
				# confront invariant mass in reco and MC - diagonal elements of confusion matrix
				.Define("MC_rho_m","Ztautau::get_MCmeson_mass(myEvent,4,true,4,true)")
				.Define("MC_a1_m","Ztautau::get_MCmeson_mass(myEvent,5,true,5,true)")
				.Define("reco_rho_m","Ztautau::get_invariant_mass(myEvent,4,true,4,true)")
				.Define("reco_a1_m","Ztautau::get_invariant_mass(myEvent,5,true,5,true)")
				.Define("rho_pull","MC_rho_m - reco_rho_m")
				.Define("a1_pull","MC_a1_m - reco_a1_m")
				  
                )
		
        #df2.Display(["MC_event","MC_event_type","weights_plus","weights_minus","found"],20).Print()
        #df2.Display(["ElAsPi_e","MuAsPi_e","RhoAsPi_e","A1AsPi_e","PiAsPi_e"],20).Print()
        return df2
       

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
        	"MC_event",
        	"event_type_reco",
        	"debug",
        	"debug_mass",
        	"pi_mask",
        	"weight_mask",
        	"pi_sgn",
			"w_plus",
			"w_minus",
			"mu_sgn",
			"el_sgn",
            "RP_thrustcostheta",
			"RP_thrustphi",
        	"ElAsPi_e",
        	"MuAsPi_e",
        	"RhoAsPi_e",
            "A1AsPi_e",
        	"PiAsPi_e",
        	"OtAsPi_e",
			"reco_rho_m",
			"reco_a1_m",
			"MC_rho_m",
			"MC_a1_m",
			"rho_pull",
			"a1_pull"
		
        	
        	]
        return branchList

