#Optional: List of analysis packages to load in runtime
#analysesList = ['Ztautau']

#List of processes
processList = {
    'p8_ee_Ztautau_ecm91':{},
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
#prodTag     = "FCCee/winter2023/IDEA/"
# Link to the dictonary that contains all the cross section informations etc... (mandatory)
procDict = "FCCee_procDict_winter2023_IDEA.json"

#output directory
outputDir = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/angular_analysis/"
#input directory
inputDir    = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"

#Optional
#nCPUS       = 1
#runBatch    = True
#batchQueue = "espresso"
#compGroup = "group_u_FCC.local_gen"

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
				
				.Filter("abs(RP_thrustcostheta0)<0.95")
                
				
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
				.Define("Pion0", "Ztautau::sel_pions_id(rps,1)")
				
				#####
				# EVENTS IDENTIFICATION
				#####
				.Define("myEvent","Ztautau::myget_event(Muon0,Electron0,Pion0,Photon0,rps,RP_thrustangle,Particle,Particle1,rp2mc_idx,EVT_thrust)")


				# Create 20 pre-binned columns for angular analysis (costheta from -0.95 to +0.95 at step 0.095)
				.Define("pi_sgn_bin0","Ztautau::get_reco_x_theta(myEvent,3,false,3,true,-0.95,-0.855)")
				.Define("pi_sgn_bin1","Ztautau::get_reco_x_theta(myEvent,3,false,3,true,-0.855,-0.760)")
				.Define("pi_sgn_bin2","Ztautau::get_reco_x_theta(myEvent,3,false,3,true,-0.760,-0.665)")
				.Define("pi_sgn_bin3","Ztautau::get_reco_x_theta(myEvent,3,false,3,true,-0.665,-0.570)")
				.Define("pi_sgn_bin4","Ztautau::get_reco_x_theta(myEvent,3,false,3,true,-0.570,-0.475)")
				.Define("pi_sgn_bin5","Ztautau::get_reco_x_theta(myEvent,3,false,3,true,-0.475,-0.380)")
				.Define("pi_sgn_bin6","Ztautau::get_reco_x_theta(myEvent,3,false,3,true,-0.380,-0.285)")
				.Define("pi_sgn_bin7","Ztautau::get_reco_x_theta(myEvent,3,false,3,true,-0.285,-0.190)")
				.Define("pi_sgn_bin8","Ztautau::get_reco_x_theta(myEvent,3,false,3,true,-0.190,-0.095)")
				.Define("pi_sgn_bin9","Ztautau::get_reco_x_theta(myEvent,3,false,3,true,-0.095,0.000)")
				.Define("pi_sgn_bin10","Ztautau::get_reco_x_theta(myEvent,3,false,3,true,0.000,0.095)")
				.Define("pi_sgn_bin11","Ztautau::get_reco_x_theta(myEvent,3,false,3,true,0.095,0.190)")
				.Define("pi_sgn_bin12","Ztautau::get_reco_x_theta(myEvent,3,false,3,true,0.190,0.285)")
				.Define("pi_sgn_bin13","Ztautau::get_reco_x_theta(myEvent,3,false,3,true,0.285,0.380)")
				.Define("pi_sgn_bin14","Ztautau::get_reco_x_theta(myEvent,3,false,3,true,0.380,0.475)")
				.Define("pi_sgn_bin15","Ztautau::get_reco_x_theta(myEvent,3,false,3,true,0.475,0.570)")
				.Define("pi_sgn_bin16","Ztautau::get_reco_x_theta(myEvent,3,false,3,true,0.570,0.665)")
				.Define("pi_sgn_bin17","Ztautau::get_reco_x_theta(myEvent,3,false,3,true,0.665,0.760)")
				.Define("pi_sgn_bin18","Ztautau::get_reco_x_theta(myEvent,3,false,3,true,0.760,0.855)")
				.Define("pi_sgn_bin19","Ztautau::get_reco_x_theta(myEvent,3,false,3,true,0.855,0.950)")

				.Define("el_sgn_bin0","Ztautau::get_reco_x_theta(myEvent,2,false,2,true,-0.95,-0.855)")
				.Define("el_sgn_bin1","Ztautau::get_reco_x_theta(myEvent,2,false,2,true,-0.855,-0.760)")
				.Define("el_sgn_bin2","Ztautau::get_reco_x_theta(myEvent,2,false,2,true,-0.760,-0.665)")
				.Define("el_sgn_bin3","Ztautau::get_reco_x_theta(myEvent,2,false,2,true,-0.665,-0.570)")
				.Define("el_sgn_bin4","Ztautau::get_reco_x_theta(myEvent,2,false,2,true,-0.570,-0.475)")
				.Define("el_sgn_bin5","Ztautau::get_reco_x_theta(myEvent,2,false,2,true,-0.475,-0.380)")
				.Define("el_sgn_bin6","Ztautau::get_reco_x_theta(myEvent,2,false,2,true,-0.380,-0.285)")
				.Define("el_sgn_bin7","Ztautau::get_reco_x_theta(myEvent,2,false,2,true,-0.285,-0.190)")
				.Define("el_sgn_bin8","Ztautau::get_reco_x_theta(myEvent,2,false,2,true,-0.190,-0.095)")
				.Define("el_sgn_bin9","Ztautau::get_reco_x_theta(myEvent,2,false,2,true,-0.095,0.000)")
				.Define("el_sgn_bin10","Ztautau::get_reco_x_theta(myEvent,2,false,2,true,0.000,0.095)")
				.Define("el_sgn_bin11","Ztautau::get_reco_x_theta(myEvent,2,false,2,true,0.095,0.190)")
				.Define("el_sgn_bin12","Ztautau::get_reco_x_theta(myEvent,2,false,2,true,0.190,0.285)")
				.Define("el_sgn_bin13","Ztautau::get_reco_x_theta(myEvent,2,false,2,true,0.285,0.380)")
				.Define("el_sgn_bin14","Ztautau::get_reco_x_theta(myEvent,2,false,2,true,0.380,0.475)")
				.Define("el_sgn_bin15","Ztautau::get_reco_x_theta(myEvent,2,false,2,true,0.475,0.570)")
				.Define("el_sgn_bin16","Ztautau::get_reco_x_theta(myEvent,2,false,2,true,0.570,0.665)")
				.Define("el_sgn_bin17","Ztautau::get_reco_x_theta(myEvent,2,false,2,true,0.665,0.760)")
				.Define("el_sgn_bin18","Ztautau::get_reco_x_theta(myEvent,2,false,2,true,0.760,0.855)")
				.Define("el_sgn_bin19","Ztautau::get_reco_x_theta(myEvent,2,false,2,true,0.855,0.950)")

				.Define("mu_sgn_bin0","Ztautau::get_reco_x_theta(myEvent,1,false,1,true,-0.95,-0.855)")
				.Define("mu_sgn_bin1","Ztautau::get_reco_x_theta(myEvent,1,false,1,true,-0.855,-0.760)")
				.Define("mu_sgn_bin2","Ztautau::get_reco_x_theta(myEvent,1,false,1,true,-0.760,-0.665)")
				.Define("mu_sgn_bin3","Ztautau::get_reco_x_theta(myEvent,1,false,1,true,-0.665,-0.570)")
				.Define("mu_sgn_bin4","Ztautau::get_reco_x_theta(myEvent,1,false,1,true,-0.570,-0.475)")
				.Define("mu_sgn_bin5","Ztautau::get_reco_x_theta(myEvent,1,false,1,true,-0.475,-0.380)")
				.Define("mu_sgn_bin6","Ztautau::get_reco_x_theta(myEvent,1,false,1,true,-0.380,-0.285)")
				.Define("mu_sgn_bin7","Ztautau::get_reco_x_theta(myEvent,1,false,1,true,-0.285,-0.190)")
				.Define("mu_sgn_bin8","Ztautau::get_reco_x_theta(myEvent,1,false,1,true,-0.190,-0.095)")
				.Define("mu_sgn_bin9","Ztautau::get_reco_x_theta(myEvent,1,false,1,true,-0.095,0.000)")
				.Define("mu_sgn_bin10","Ztautau::get_reco_x_theta(myEvent,1,false,1,true,0.000,0.095)")
				.Define("mu_sgn_bin11","Ztautau::get_reco_x_theta(myEvent,1,false,1,true,0.095,0.190)")
				.Define("mu_sgn_bin12","Ztautau::get_reco_x_theta(myEvent,1,false,1,true,0.190,0.285)")
				.Define("mu_sgn_bin13","Ztautau::get_reco_x_theta(myEvent,1,false,1,true,0.285,0.380)")
				.Define("mu_sgn_bin14","Ztautau::get_reco_x_theta(myEvent,1,false,1,true,0.380,0.475)")
				.Define("mu_sgn_bin15","Ztautau::get_reco_x_theta(myEvent,1,false,1,true,0.475,0.570)")
				.Define("mu_sgn_bin16","Ztautau::get_reco_x_theta(myEvent,1,false,1,true,0.570,0.665)")
				.Define("mu_sgn_bin17","Ztautau::get_reco_x_theta(myEvent,1,false,1,true,0.665,0.760)")
				.Define("mu_sgn_bin18","Ztautau::get_reco_x_theta(myEvent,1,false,1,true,0.760,0.855)")
				.Define("mu_sgn_bin19","Ztautau::get_reco_x_theta(myEvent,1,false,1,true,0.855,0.950)")

				.Define("rho_sgn_bin0","Ztautau::get_reco_x_theta(myEvent,4,false,4,true,-0.95,-0.855)")
				.Define("rho_sgn_bin1","Ztautau::get_reco_x_theta(myEvent,4,false,4,true,-0.855,-0.760)")
				.Define("rho_sgn_bin2","Ztautau::get_reco_x_theta(myEvent,4,false,4,true,-0.760,-0.665)")
				.Define("rho_sgn_bin3","Ztautau::get_reco_x_theta(myEvent,4,false,4,true,-0.665,-0.570)")
				.Define("rho_sgn_bin4","Ztautau::get_reco_x_theta(myEvent,4,false,4,true,-0.570,-0.475)")
				.Define("rho_sgn_bin5","Ztautau::get_reco_x_theta(myEvent,4,false,4,true,-0.475,-0.380)")
				.Define("rho_sgn_bin6","Ztautau::get_reco_x_theta(myEvent,4,false,4,true,-0.380,-0.285)")
				.Define("rho_sgn_bin7","Ztautau::get_reco_x_theta(myEvent,4,false,4,true,-0.285,-0.190)")
				.Define("rho_sgn_bin8","Ztautau::get_reco_x_theta(myEvent,4,false,4,true,-0.190,-0.095)")
				.Define("rho_sgn_bin9","Ztautau::get_reco_x_theta(myEvent,4,false,4,true,-0.095,0.000)")
				.Define("rho_sgn_bin10","Ztautau::get_reco_x_theta(myEvent,4,false,4,true,0.000,0.095)")
				.Define("rho_sgn_bin11","Ztautau::get_reco_x_theta(myEvent,4,false,4,true,0.095,0.190)")
				.Define("rho_sgn_bin12","Ztautau::get_reco_x_theta(myEvent,4,false,4,true,0.190,0.285)")
				.Define("rho_sgn_bin13","Ztautau::get_reco_x_theta(myEvent,4,false,4,true,0.285,0.380)")
				.Define("rho_sgn_bin14","Ztautau::get_reco_x_theta(myEvent,4,false,4,true,0.380,0.475)")
				.Define("rho_sgn_bin15","Ztautau::get_reco_x_theta(myEvent,4,false,4,true,0.475,0.570)")
				.Define("rho_sgn_bin16","Ztautau::get_reco_x_theta(myEvent,4,false,4,true,0.570,0.665)")
				.Define("rho_sgn_bin17","Ztautau::get_reco_x_theta(myEvent,4,false,4,true,0.665,0.760)")
				.Define("rho_sgn_bin18","Ztautau::get_reco_x_theta(myEvent,4,false,4,true,0.760,0.855)")
				.Define("rho_sgn_bin19","Ztautau::get_reco_x_theta(myEvent,4,false,4,true,0.855,0.950)")

				  
                )
		
        return df2
       

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
			"pi_sgn_bin0", "pi_sgn_bin1", "pi_sgn_bin2", "pi_sgn_bin3", "pi_sgn_bin4",
			"pi_sgn_bin5", "pi_sgn_bin6", "pi_sgn_bin7", "pi_sgn_bin8", "pi_sgn_bin9",
			"pi_sgn_bin10", "pi_sgn_bin11", "pi_sgn_bin12", "pi_sgn_bin13", "pi_sgn_bin14",
			"pi_sgn_bin15", "pi_sgn_bin16", "pi_sgn_bin17", "pi_sgn_bin18", "pi_sgn_bin19",
			"el_sgn_bin0", "el_sgn_bin1", "el_sgn_bin2", "el_sgn_bin3", "el_sgn_bin4",
			"el_sgn_bin5", "el_sgn_bin6", "el_sgn_bin7", "el_sgn_bin8", "el_sgn_bin9",
			"el_sgn_bin10", "el_sgn_bin11", "el_sgn_bin12", "el_sgn_bin13", "el_sgn_bin14",
			"el_sgn_bin15", "el_sgn_bin16", "el_sgn_bin17", "el_sgn_bin18", "el_sgn_bin19",
			"mu_sgn_bin0", "mu_sgn_bin1", "mu_sgn_bin2", "mu_sgn_bin3", "mu_sgn_bin4",
			"mu_sgn_bin5", "mu_sgn_bin6", "mu_sgn_bin7", "mu_sgn_bin8", "mu_sgn_bin9",
			"mu_sgn_bin10", "mu_sgn_bin11", "mu_sgn_bin12", "mu_sgn_bin13", "mu_sgn_bin14",
			"mu_sgn_bin15", "mu_sgn_bin16", "mu_sgn_bin17", "mu_sgn_bin18", "mu_sgn_bin19",
			"rho_sgn_bin0", "rho_sgn_bin1", "rho_sgn_bin2", "rho_sgn_bin3", "rho_sgn_bin4",
			"rho_sgn_bin5", "rho_sgn_bin6", "rho_sgn_bin7", "rho_sgn_bin8", "rho_sgn_bin9",
			"rho_sgn_bin10", "rho_sgn_bin11", "rho_sgn_bin12", "rho_sgn_bin13", "rho_sgn_bin14",
			"rho_sgn_bin15", "rho_sgn_bin16", "rho_sgn_bin17", "rho_sgn_bin18", "rho_sgn_bin19"
            
		]
        return branchList

