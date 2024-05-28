import os, sys
import subprocess
import ROOT
import argparse
import time
from utils.CMShelper import *

parser = argparse.ArgumentParser(description='Skim full tuple.')
parser.add_argument('-i', '--input', required=False, type=str, default="/VBFHHto2B2Tau_CV-1_C2V-1_C3-1_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM", help="output directory name")
parser.add_argument('-o', '--output', required=False, type=str, default="./OutPut", help="output directory name")
args = parser.parse_args()

ROOT.gInterpreter.Declare('#include "./interface/HHbbtautau_PreSkim_jin_reco.h"')
ROOT.EnableImplicitMT(10)

os.makedirs(args.output, exist_ok=True)
data_base = dasgoclient(args.input)
#data_frame = remote_reading([data_base[0]]) 
data_frame = remote_reading(data_base) 
# Licheng: the problem left here is to separate database into batches when facing a large dataset

print( "We proceed {0} events in this batch\n".format( str(data_frame.Count().GetValue()) ) )

data_frame = data_frame.Filter("nFatJet >1")
print( "{0} events left after [nFatJet >1] in this batch.\n".format( str(data_frame.Count().GetValue()) ) )

data_frame = data_frame.Filter("min_result(FatJet_pt) > 180").Filter("max_result(FatJet_eta) < 2.4").Filter("min_result(FatJet_eta) > -2.4")
print( "{0} events left after [FatJet_pt > 180 && |FatJet_eta| < 2.4] in this batch.\n".format( str(data_frame.Count().GetValue()) ) )

data_frame = data_frame.Filter("veto_e_muon(Electron_pt, Electron_eta, Electron_dxy, Electron_dz, 20, 2.5, 0.05, 0.2) == 0")
#data_frame = data_frame.Filter("(max_result(Electron_pt) < 10)||(min_result(Electron_eta) > 2.4)||(max_result(Electron_eta) < -2.4)||\
#                                (min_result(Electron_dxy) > 0.05)||(max_result(Electron_dxy) < -0.05)||\
#                                (min_result(Electron_dz) > 0.2)||(max_result(Electron_dz) < -0.2)")

print( "{0} events left after [e veto] in this batch.\n".format( str(data_frame.Count().GetValue()) ) )

data_frame = data_frame.Filter("veto_e_muon(Muon_pt, Muon_eta, Muon_dxy, Muon_dz, 10, 2.4, 0.05, 0.2) == 0")
#data_frame = data_frame.Filter("(max_result(Muon_pt) < 20)||(min_result(Muon_eta) > 2.5)||(max_result(Muon_eta) < -2.5)||\
#                                (min_result(Muon_dxy) > 0.05)||(max_result(Muon_dxy) < -0.05)||\
#                                (min_result(Muon_dz) > 0.2)||(max_result(Muon_dz) < -0.2)")

print( "{0} events left after [Muon veto] in this batch.\n".format( str(data_frame.Count().GetValue()) ) )


data_frame = data_frame.Define("tautauJet_id", "find_tautauJet_from_best_PNetxtt(FatJet_particleNet_XttVsQCD)")
data_frame = data_frame.Define("mass_vistautau","FatJet_msoftdrop[tautauJet_id]*FatJet_particleNet_massCorr[tautauJet_id]")

'''data_frame = data_frame.Define("n_boostedtau", "find_2tau(FatJet_eta, FatJet_phi, boostedTau_eta, boostedTau_phi, tautauJet_id, boostedTau_pt)[0]").Filter("n_boostedtau > 1")

print( "{0} events left after [n_boostedtau > 1] in this batch.\n".format( str(data_frame.Count().GetValue()) ) )


data_frame = data_frame.Define("tau0_id", "find_2tau(FatJet_eta, FatJet_phi, boostedTau_eta, boostedTau_phi, tautauJet_id, boostedTau_pt)[1]")
data_frame = data_frame.Define("tau1_id", "find_2tau(FatJet_eta, FatJet_phi, boostedTau_eta, boostedTau_phi, tautauJet_id, boostedTau_pt)[2]")

data_frame = data_frame.Define("mass_CA", "M_tautau_CA(boostedTau_eta, boostedTau_phi, boostedTau_mass, \
                                boostedTau_pt, tau0_id, tau1_id, MET_pt, MET_phi, mass_vistautau)")'''


data_frame = data_frame.Define("n_boostedtau", "find_2tau(FatJet_eta, FatJet_phi, Tau_eta, Tau_phi, tautauJet_id, Tau_pt)[0]").Filter("n_boostedtau > 1")

print( "{0} events left after [n_boostedtau > 1] in this batch.\n".format( str(data_frame.Count().GetValue()) ) )


data_frame = data_frame.Define("tau0_id", "find_2tau(FatJet_eta, FatJet_phi, Tau_eta, Tau_phi, tautauJet_id, Tau_pt)[1]")
data_frame = data_frame.Define("tau1_id", "find_2tau(FatJet_eta, FatJet_phi, Tau_eta, boostedTau_phi, tautauJet_id, Tau_pt)[2]")

data_frame = data_frame.Define("mass_CA", "M_tautau_CA(Tau_eta, Tau_phi, Tau_mass, \
                                Tau_pt, tau0_id, tau1_id, MET_pt, MET_phi, mass_vistautau)")

# number of gen taus(events with less than two taus or 2 visibles taus are filtered out)


print( "{0} events left in this batch.\n".format( str(data_frame.Count().GetValue()) ) )

skim_branches = ["mass_vistautau","tautauJet_id", "n_boostedtau",
                
                "boostedTau_mass", "boostedTau_phi", "boostedTau_eta", "boostedTau_pt",
                
                "FatJet_eta", "FatJet_mass", "FatJet_phi", "FatJet_pt",
                
                "tau0_id", "tau1_id", "mass_CA"]

output_file = args.output + "/" + "HHbbtautau_PreSkim_test_10w_withoutboosted.root"
data_frame.Snapshot("Events", output_file, skim_branches)

