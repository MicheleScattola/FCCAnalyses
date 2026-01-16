#Optional: List of analysis packages to load in runtime
#analysesList = ['Ztautau']

#List of processes
processList = {
    'p8_ee_Ztautau_ecm91':{'fraction':0.03},
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
nCPUS       = 1
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
				.Define("pions_charged_ids", "Ztautau::sel_pions_id(rps,1)")
				
				#####
				# EVENTS IDENTIFICATION
				#####
				.Define("myEvent","Ztautau::myget_event(Muon0,Electron0,pions_charged_ids,Photon0,rps,RP_thrustangle,Particle,Particle1)")
                
				.Define("event_type_reco","Ztautau::get_type_safe(myEvent)")
				
				#####
				# MC IDENTIFICATION
				#####
                .Define("MC_event","RVec<int> {myEvent[0].m_MCtype,myEvent[1].m_MCtype}")
				
				# parameters are: event collection, MC type, reco type, masscheck
                
				# background sources
                
				# pi confusion
				.Define("ElAsPi_e","Ztautau::get_hadron_e(myEvent,2,true,3,true)")
				.Define("MuAsPi_e","Ztautau::get_hadron_e(myEvent,1,true,3,true)")
				.Define("RhoAsPi_e","Ztautau::get_hadron_e(myEvent,4,true,3,true)")
				.Define("PiAsPi_e","Ztautau::get_hadron_e(myEvent,3,true,3,true)")
				.Define("A1AsPi_e","Ztautau::get_hadron_e(myEvent,5,true,3,true)")
                
				# rho confusion
                .Define("PiAsRho_e","Ztautau::get_hadron_e(myEvent,3,true,4,true)")
                .Define("RhoAsRho_e","Ztautau::get_hadron_e(myEvent,4,true,4,true)")
                .Define("A1AsRho_e","Ztautau::get_hadron_e(myEvent,5,true,4,true)")
                .Define("AllRho_e","Ztautau::get_hadron_e(myEvent,4,false,4,true)")
                
				  
                )
		
        return df2
       

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
        	"MC_event",
        	"event_type_reco",
        	"ElAsPi_e",
        	"MuAsPi_e",
        	"RhoAsPi_e",
            "A1AsPi_e",
        	"PiAsPi_e",
            "PiAsRho_e",
			"RhoAsRho_e",
			"A1AsRho_e",
			"AllRho_e"
            
		
        	
        	]
        return branchList

