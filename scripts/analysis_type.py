#Optional: List of analysis packages to load in runtime
#analysesList = ['Ztautau']

#List of processes
processList = {
    'p8_ee_Ztautau_ecm91':{'fraction':0.01},
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
				.Define("event_type_all","Ztautau::get_event( muons,muons_costheta, electrons, electrons_costheta, pions_charged, pi_costheta, photons, photons_costheta,5.)")
				.Define("event_type_reco","RVec<int> {event_type_all[0].type,event_type_all[1].type}")
				
				#####
				# MC IDENTIFICATION
				#####
                .Define("MC_event","Ztautau::classify_mc_event(Particle,Particle1)")
				
				# my particle collections
				.Define("myEl","Ztautau::buildRPTruthCollection(rps,Particle,Particle0, rps_MC_index,RP_thrustangle,Electron0,MC_event,event_type_reco)")
				.Define("myMu","Ztautau::buildRPTruthCollection(rps,Particle,Particle0, rps_MC_index,RP_thrustangle,Muon0,MC_event,event_type_reco)")
				.Define("myPi","Ztautau::buildRPTruthCollection(rps,Particle,Particle0, rps_MC_index,RP_thrustangle,pions_charged_ids,MC_event,event_type_reco)")
				.Define("myPh","Ztautau::buildRPTruthCollection(rps,Particle,Particle0, rps_MC_index,RP_thrustangle,Photon0,MC_event,event_type_reco)")
				
				# photon involved
				#.Define("RhoAsPi_phe","Ztautau::confusion_e(myPh,4,true,3,true,22,true,0)")
				#.Define("OtAsRho_phe","Ztautau::confusion_e(myPh,0,true,4,true,22,true,0)")
				#.Define("RhoAsOt_phe","Ztautau::confusion_e(myPh,4,true,0,true,22,true,0)")
				
				#.Define("RhoAsPi_pht","Ztautau::confusion_theta(myPh,4,true,3,true,22,true)")
				#.Define("RhoAsOt_pht","Ztautau::confusion_theta(myPh,4,true,0,true,22,true)")
				#.Define("OtAsRho_pht","Ztautau::confusion_theta(myPh,0,true,4,true,22,true)")
				
				
				# Miss ID as 1prong pions
				.Define("ElAsPi_t","Ztautau::confusion_theta(myPi,2,true,3,true,11,true)")
				.Define("MuAsPi_t","Ztautau::confusion_theta(myPi,1,true,3,true,13,true)")
				.Define("RhoAsPi_t","Ztautau::confusion_theta(myPi,4,true,3,true,111,false)")
				.Define("PiAsPi_t","Ztautau::confusion_theta(myPi,3,true,3,true,211,false)")
				.Define("A1AsPi_t","Ztautau::confusion_theta(myPi,5,true,3,true,211,false)")
                .Define("OtAsPi_t","Ztautau::confusion_theta(myPi,0,true,3,true,211,false)")
				
				.Define("ElAsPi_e","Ztautau::confusion_e(myPi,2,true,3,true,11,true,0)")
				.Define("MuAsPi_e","Ztautau::confusion_e(myPi,1,true,3,true,13,true,0)")
				.Define("RhoAsPi_e","Ztautau::confusion_e(myPi,4,true,3,true,111,false,0)")
				.Define("PiAsPi_e","Ztautau::confusion_e(myPi,3,true,3,true,211,false,0)")
				.Define("OtAsPi_e","Ztautau::confusion_e(myPi,0,true,3,true,211,false,0)")
                .Define("A1AsPi_e","Ztautau::confusion_e(myPi,5,true,3,true,211,false,0)")
				
				  
                )
		
        #df2.Display(["myEl_e","myMu_e","MC_event","event_type_reco"],20).Print()
        
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
        	"ElAsPi_t",
        	"MuAsPi_t",
        	"RhoAsPi_t",
            "A1AsPi_e",
            "A1AsPi_t",
        	"PiAsPi_e",
        	"PiAsPi_t",
        	"OtAsPi_e",
        	"OtAsPi_t"
        	
        	]
        return branchList

