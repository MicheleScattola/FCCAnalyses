#Optional: List of analysis packages to load in runtime
#analysesList = ['Ztautau']

#List of processes
processList = {
    'p8_ee_Ztautau_ecm91':{'fraction':0.02,'chunks':4},
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
#prodTag     = "FCCee/winter2023/IDEA/"
# Link to the dictonary that contains all the cross section informations etc... (mandatory)
procDict = "FCCee_procDict_winter2023_IDEA.json"

#output directory
outputDir = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/templates/"
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
				.Define("rp2mc_idx",
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
                
				# pi signal and weights
				.Define("pi_sgn","Ztautau::get_hadron_e(myEvent,3,false,3,true)/45.5")
				.Define("w_plus_pi","Ztautau::get_weights(1,myEvent,3,false,3,true,true,false)")
				.Define("w_minus_pi","Ztautau::get_weights(-1,myEvent,3,false,3,true,true,false)")
				# rho signal
                .Define("rho_sgn","Ztautau::get_omega_rho(myEvent,4,false,4,true,true,false)")
                .Define("w_plus_rho","Ztautau::get_weights(1,myEvent,4,false,4,true,true,false)")
				.Define("w_minus_rho","Ztautau::get_weights(-1,myEvent,4,false,4,true,true,false)")
                
				# lep signal and weights
                .Define("el_sgn","Ztautau::get_lepton_e(myEvent,2,false,2,true,true,false)/45.5")
                .Define("mu_sgn","Ztautau::get_lepton_e(myEvent,1,false,1,true,true,false)/45.5")
                .Define("w_plus_el","Ztautau::get_weights(1,myEvent,2,false,2,true,true,false)")
				.Define("w_minus_el","Ztautau::get_weights(-1,myEvent,2,false,2,true,true,false)")
                .Define("w_plus_mu","Ztautau::get_weights(1,myEvent,1,false,1,true,true,false)")
				.Define("w_minus_mu","Ztautau::get_weights(-1,myEvent,1,false,1,true,true,false)")
				
				  
                )
		
        #df2.Display(["MC_event","MC_event_type","weights_plus","weights_minus","found"],20).Print()
        #df2.Display(["ElAsPi_e","MuAsPi_e","RhoAsPi_e","A1AsPi_e","PiAsPi_e"],20).Print()
        return df2
       

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
        	"pi_sgn",
			"w_plus_pi",
			"w_minus_pi",
            "el_sgn",
            "w_plus_el",
            "w_minus_el",
            "mu_sgn",
            "w_plus_mu",
            "w_minus_mu",
            "rho_sgn",
			"w_plus_rho",
            "w_minus_rho",
        	
        	]
        return branchList

