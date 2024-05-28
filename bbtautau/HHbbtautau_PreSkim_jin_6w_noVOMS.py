import os, sys
import subprocess
import ROOT
import argparse
import time
from utils.CMShelper import *

parser = argparse.ArgumentParser(description='Skim full tuple.')

parser.add_argument('-m', '--mode', required=False, type=str, default="n", help="output directory name")
parser.add_argument('-l', '--localinput', required=False, type=str, default="/eos/user/j/jinwa/bbtautau/HEPStuff/bbtautau/input/VBFHHbbtautau/294beabb-42f5-4d00-898c-9ed7e0e90049.root", help="output directory name")

parser.add_argument('-i', '--input', required=False, type=str, default="/VBFHHto2B2Tau_CV-1_C2V-1_C3-1_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM", help="output directory name")
parser.add_argument('-o', '--output', required=False, type=str, default="./OutPut", help="output directory name")

parser.add_argument('-c', '--emuoncut', required=False, type=str, default="y", help="Lichenghaoshuai")
parser.add_argument('-b', '--usingboostedtau', required=False, type=str, default="y", help="Lichenghaoshuai")


args = parser.parse_args()

ROOT.gInterpreter.Declare('#include "./interface/HHbbtautau_PreSkim_jin_reco.h"')
ROOT.EnableImplicitMT(10)

os.makedirs(args.output, exist_ok=True)

if(args.mode == "n"):
    data_frame = ROOT.RDataFrame('Events', [args.localinput])
    if("VBFHHbbtautau" in args.localinput):
        channel = "VBFHHbbtautau"
        data_frame = data_frame.Define("Weight_check", "genWeight * 1.0 / 1.0")

    elif("DYto2L-2Jets_MLL-50" in args.localinput):
        channel = "DYto2L-2Jets_MLL-50"
        data_frame = data_frame.Define("Weight_check", "genWeight * 65997137.0 / 439394.0")
    elif("DYto2L-2Jets_MLL-10to50" in args.localinput):
        channel = "DYto2L-2Jets_MLL-10to50"
        data_frame = data_frame.Define("Weight_check", "genWeight * 67681922.0 / 973875.0")
    elif("DYto2L-4Jets_MLL-10to50" in args.localinput):
        channel = "DYto2L-2Jets_MLL-50"
        data_frame = data_frame.Define("Weight_check", "genWeight * 154413937.0 / 1442259.0")
    
    elif("TTto4Q" in args.localinput):
        channel = "TTto4Q"
        data_frame = data_frame.Define("Weight_check", "genWeight * 53605620.0 / 594428.0 ")
    elif("TTto2L2Nu" in args.localinput):
        channel = "TTto2L2Nu"
        data_frame = data_frame.Define("Weight_check", "genWeight * 23778148.0 / 587160.0 ")
    elif("TTtoLNu2Q" in args.localinput):
        channel = "TTtoLNu2Q"
        data_frame = data_frame.Define("Weight_check", "genWeight * 87986940.0 / 542100.0 ")
else:
    data_base = dasgoclient(args.input)
    data_frame = remote_reading(data_base)
    channel = "DAS" 


print( "We proceed {0} events in this batch\n".format( str(data_frame.Count().GetValue()) ) )
events_proceed = str(data_frame.Count().GetValue())

data_frame = data_frame.Filter("nFatJet >1")
print( "{0} events left after [nFatJet >1] in this batch.\n".format( str(data_frame.Count().GetValue()) ) )

data_frame = data_frame.Filter("min_result(FatJet_pt) > 180").Filter("max_result(FatJet_eta) < 2.4").Filter("min_result(FatJet_eta) > -2.4")
print( "{0} events left after [FatJet_pt > 180 && |FatJet_eta| < 2.4] in this batch.\n".format( str(data_frame.Count().GetValue()) ) )

if(args.emuoncut == "y"):
    
    data_frame = data_frame.Filter("veto_e_muon(Electron_pt, Electron_eta, Electron_dxy, Electron_dz, 20, 2.5, 0.05, 0.2) == 0")
    print( "{0} events left after [e veto] in this batch.\n".format( str(data_frame.Count().GetValue()) ) )
    
    data_frame = data_frame.Filter("veto_e_muon(Muon_pt, Muon_eta, Muon_dxy, Muon_dz, 10, 2.4, 0.05, 0.2) == 0")
    print( "{0} events left after [Muon veto] in this batch.\n".format( str(data_frame.Count().GetValue()) ) )
else:
    pass

data_frame = data_frame.Define("tautauJet_id", "find_tautauJet_from_best_PNetxtt(FatJet_particleNet_XttVsQCD)")
data_frame = data_frame.Define("mass_vistautau","FatJet_msoftdrop[tautauJet_id]*FatJet_particleNet_massCorr[tautauJet_id]")
data_frame = data_frame.Define("mass_vistautau_nocorr","FatJet_msoftdrop[tautauJet_id]")

data_frame = data_frame.Define("TauTauJet_particleNet_XttVsQCD", "FatJet_particleNet_XttVsQCD[tautauJet_id]")


if(args.usingboostedtau == "y"):
    data_frame = data_frame.Define("boostedtau_vector", "find_2tau(FatJet_eta, FatJet_phi, boostedTau_eta, boostedTau_phi, tautauJet_id, boostedTau_pt)")

    data_frame = data_frame.Define("n_boostedtau", "boostedtau_vector[0]").Filter("n_boostedtau > 1")
    
    print( "{0} events left after [n_boostedtau > 1] in this batch.\n".format( str(data_frame.Count().GetValue()) ) )
    
    data_frame = data_frame.Define("tau0_id", "boostedtau_vector[1]")
    data_frame = data_frame.Define("tau1_id", "boostedtau_vector[2]")
    
    data_frame = data_frame.Define("mass_CA", "M_tautau_CA(boostedTau_eta, boostedTau_phi, boostedTau_mass, \
                                    boostedTau_pt, tau0_id, tau1_id, MET_pt, MET_phi, mass_vistautau)")

else:
    data_frame = data_frame.Define("tau_vector", "find_2tau(FatJet_eta, FatJet_phi, Tau_eta, Tau_phi, tautauJet_id, Tau_pt)")
    
    data_frame = data_frame.Define("n_boostedtau", "tau_vector[0]").Filter("n_boostedtau > 1")
    
    print( "{0} events left after [n_boostedtau > 1] in this batch.\n".format( str(data_frame.Count().GetValue()) ) )
    
    data_frame = data_frame.Define("tau0_id", "tau_vector[1]")
    data_frame = data_frame.Define("tau1_id", "tau_vector[2]")
    
    data_frame = data_frame.Define("mass_CA", "M_tautau_CA(Tau_eta, Tau_phi, Tau_mass, \
                                    Tau_pt, tau0_id, tau1_id, MET_pt, MET_phi, mass_vistautau)")

# number of gen taus(events with less than two taus or 2 visibles taus are filtered out)


print( "{0} events left in this batch.\n".format( str(data_frame.Count().GetValue()) ) )

skim_branches = ["mass_vistautau", "mass_vistautau_nocorr", "mass_CA", "genWeight", "Weight_check",

                "tautauJet_id", "TauTauJet_particleNet_XttVsQCD", "n_boostedtau",
                
                "boostedTau_mass", "boostedTau_phi", "boostedTau_eta", "boostedTau_pt",
                
                "FatJet_eta", "FatJet_mass", "FatJet_phi", "FatJet_pt",
                
                "tau0_id", "tau1_id"]


output_file = args.output + "/" + channel + "_" + events_proceed + "_usingboostedtau_" + args.usingboostedtau + "_emcut_" + args.emuoncut + ".root"

data_frame.Snapshot("Events", output_file, skim_branches)

print( "Done")
