#Optional: List of analysis packages to load in runtime
#analysesList = ['Ztautau']

#List of processes
processList = {
    'p8_ee_Ztautau_ecm91':{'fraction':0.1 , 'chunks':20},
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
#prodTag     = "FCCee/winter2023/IDEA/"
# Link to the dictonary that contains all the cross section informations etc... (mandatory)
procDict = "FCCee_procDict_winter2023_IDEA.json"

#output directory
outputDir = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/"
outputName = "p8_ee_Ztautau_ecm91"
#input directory
inputDir    = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"

#Optional
nCPUS       = 1
runBatch    = True
batchQueue = "espresso"
compGroup = "group_u_FCC.local_gen"

#Optional test file , run with --test
testFile = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/p8_ee_Ztautau_ecm91/events_080694422.root"

#operations on the TTreout.push_back(e.m_type);out.push_back(e.m_type);e
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
                .Define("rp2mc_idx","ReconstructedParticle2MC::getRP2MC_index(MCRecoAssociations0, MCRecoAssociations1, rps)")
                
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
				
				.Filter("abs(RP_thrustcostheta)<0.95")
				
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

				# Electrons
				.Alias("Electron0", "Electron#0.index")

				# Photons
				.Alias("Photon0", "Photon#0.index")
				
				# defining pions as charged hadrons with mass selection
				# selecting candidates (possibly mistaken with a K+ )
				.Define("Pion0", "Ztautau::sel_pions_id(rps,1)")
				
				#####
				# EVENTS IDENTIFICATION
				#####
				.Define("myEvent","Ztautau::myget_event(Muon0,Electron0,Pion0,Photon0,rps,RP_thrustangle,Particle,Particle1,rp2mc_idx,RP_thrustcostheta,RP_thrustphi)")
                
				.Define("event_type_reco","Ztautau::get_type_safe(myEvent,true)")
                .Define("event_type_reco_noMassCheck","Ztautau::get_type_safe(myEvent,false)")
				
				#####
				# MC IDENTIFICATION
				#####
                .Define("MC_event","RVec<int> {myEvent[0].mc_type,myEvent[1].mc_type}")
                
				# masks for later cuts and selections
				#.Define("pi_mask", "Ztautau::get_pi_mask(myEvent)")
				#.Define("weight_mask", "Ztautau::get_weight_mask(myEvent)")
                
				#APPLYING INV MASS CHECK:
				# pi signal
				.Define("pi_sgn","Ztautau::get_hadron_e(myEvent,3,false,3,true)/45.5")
				# lepton signals
                .Define("el_sgn","Ztautau::get_lepton_e(myEvent,2,false,2,true,true,false)/45.5")
                .Define("mu_sgn","Ztautau::get_lepton_e(myEvent,1,false,1,true,true,false)/45.5")
				# rho signal
                .Define("rho_sgn","Ztautau::get_omega_rho(myEvent,4,false,4,true,true,false)")

				# MC invariant mass
				.Define("MC_rho_m","Ztautau::get_MCdaughter_mass(myEvent,4,true,4,false)")
				.Define("MC_a1_m","Ztautau::get_MCdaughter_mass(myEvent,5,true,5,false)")
                # Reco invariant mass
				.Define("reco_rho_m","Ztautau::get_invariant_mass(myEvent,4,false,4,true)")
				.Define("reco_a1_m","Ztautau::get_invariant_mass(myEvent,5,false,5,true)")
				.Define("reco_pi_m","Ztautau::get_invariant_mass(myEvent,3,false,3,true)")
                # confront invariant mass in reco and MC - diagonal elements of confusion matrix
				.Define("rho_pull","Ztautau::get_mass_pull(myEvent,4,true,4,true)")
				.Define("a1_pull","Ztautau::get_mass_pull(myEvent,5,true,5,true)")
                
				# debug mass 
                .Define("rho_reject","Ztautau::get_debug_mass(myEvent,4)")
                .Define("a1_reject","Ztautau::get_debug_mass(myEvent,5)")

				  
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
            "event_type_reco_noMassCheck",
        	"pi_sgn",
			"mu_sgn",
			"el_sgn",
            "rho_sgn",
            "RP_thrustcostheta",
			"RP_thrustphi",
			"reco_rho_m",
			"reco_a1_m",
			"reco_pi_m",
			"MC_rho_m",
			"MC_a1_m",
			"rho_pull",
			"a1_pull",
            "rho_reject",
			"a1_reject"
		
        	
        	]
        return branchList

