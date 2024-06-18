import os, sys
import subprocess
import ROOT
import argparse
import time
from utils.CMShelper import *

parser = argparse.ArgumentParser(description='Skim full tuple.')

parser.add_argument('-m', '--mode', required=False, type=str, default="n", help="output directory name")
parser.add_argument('-l', '--localinput', required=False, type=str, default="/eos/user/j/jinwa/bbtautau/HEPStuff/bbtautau/input/VBFbbttCV1p21/61300cf4-2944-491e-adf1-309fb1f2866e.root", help="output directory name")

parser.add_argument('-i', '--input', required=False, type=str, default="/VBFHHto2B2Tau_CV-1_C2V-1_C3-1_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM", help="output directory name")
parser.add_argument('-o', '--output', required=False, type=str, default="./OutPut", help="output directory name")

parser.add_argument('-c', '--emuoncut', required=False, type=str, default="y", help="Lichenghaoshuai")
parser.add_argument('-b', '--usingboostedtau', required=False, type=str, default="n", help="Lichenghaoshuai")

parser.add_argument('-z', '--zero_mass', required=False, type=str, default="-1", help="Lichenghaoshuai")
parser.add_argument('-g', '--gen_tau_match', required=False, type=str, default="y", help="Lichenghaoshuai")
parser.add_argument('-r', '--reco_level_Ana', required=False, type=str, default="y", help="lchs")
parser.add_argument('-f', '--FJ_PT_cut', required=False, type=str, default="180", help="lchs")

args = parser.parse_args()

ROOT.gInterpreter.Declare('#include "./interface/HHbbtautau_PreSkim_jin_efficiency_check.h"')
ROOT.EnableImplicitMT(10)

os.makedirs(args.output, exist_ok=True)

if(args.mode == "n"):
    data_frame = ROOT.RDataFrame('Events', [args.localinput])
    if("VBFHHbbtautau" in args.localinput):
        channel = "VBFHHbbtautau"
        data_frame = data_frame.Define("Weight_check", "genWeight * 1.0 / 1.0")
    elif("VBFbbttCV1p21" in args.localinput):
        channel = "VBFbbttCV1p21"
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


if(args.reco_level_Ana == "n"):

    print("Selecting Gen Taus...\n")
    data_frame = data_frame.Filter("nGenVisTau == 2")
    # number of gen taus(events with less than two taus or 2 visibles taus are filtered out)
    print( "{0} events left in this batch after 2 Gen Visible Tau selection.\n".format( str(data_frame.Count().GetValue()) ) )

    data_frame = data_frame.Filter("max_result(FatJet_pt) > "+args.FJ_PT_cut)
    print( "{0} events left after".format( str(data_frame.Count().GetValue()) ) + "FatJet_pt > {0} in this batch.\n".format(args.FJ_PT_cut ) )
    
    data_frame = data_frame.Filter("nFatJet >1")
    print( "{0} events left after [nFatJet >1] in this batch.\n".format( str(data_frame.Count().GetValue()) ) )
    


    data_frame = data_frame.Define("Matched_FatJet","matching_FatJet(FatJet_eta, FatJet_phi, GenVisTau_eta, GenVisTau_phi)")

    data_frame = data_frame.Define("N_Matched_FatJet","Matched_FatJet[0]")
    data_frame = data_frame.Define("Matched_FatJet_id","Matched_FatJet[1]")
    data_frame = data_frame.Define("N_Matched_1Vistau_FetJet","matching_FatJet_number(FatJet_eta, FatJet_phi, GenVisTau_eta, GenVisTau_phi)[1]")
    if(args.gen_tau_match == "y"):
        data_frame = data_frame.Filter("N_Matched_FatJet > 0")
        print( "{0} events left after [1 FatJet matched with 2 GenVisTau] in this batch.\n".format( str(data_frame.Count().GetValue()) ) )
    elif(args.gen_tau_match == "n"):
        pass

    data_frame = data_frame.Define("PNet_ttvsQCD_Matched","find_tautauJet_from_best_PNetxtt(FatJet_particleNet_XttVsQCD) == Matched_FatJet[1]")


elif(args.reco_level_Ana == "y"):
    
    data_frame = data_frame.Filter("nFatJet >1")
    print( "{0} events left after [nFatJet >1] in this batch.\n".format( str(data_frame.Count().GetValue()) ) )

    data_frame = data_frame.Filter("max_result(FatJet_pt) > 180")
    print( "{0} events left after FatJet_pt > 180  in this batch.\n".format( str(data_frame.Count().GetValue()) ) )

    if(args.emuoncut == "y"):
        
        data_frame = data_frame.Filter("veto_e_muon(Electron_pt, Electron_eta, Electron_dxy, Electron_dz, 20, 2.5, 0.05, 0.2) == 0")
        print( "{0} events left after [e veto] in this batch.\n".format( str(data_frame.Count().GetValue()) ) )
        
        data_frame = data_frame.Filter("veto_e_muon(Muon_pt, Muon_eta, Muon_dxy, Muon_dz, 10, 2.4, 0.05, 0.2) == 0")
        print( "{0} events left after [Muon veto] in this batch.\n".format( str(data_frame.Count().GetValue()) ) )
    else:
        pass

    data_frame = data_frame.Define("Matched_FatJet_id", "find_tautauJet_from_best_PNetxtt(FatJet_particleNet_XttVsQCD)")
    
    data_frame = data_frame.Define("GenTau_vector", "find_2tau(FatJet_eta, FatJet_phi, GenVisTau_eta, GenVisTau_phi, Matched_FatJet_id, GenVisTau_pt)")
    data_frame = data_frame.Define("N_GenTau_Matched", "GenTau_vector[0]")


data_frame = data_frame.Define("mass_tauBoostedJet","FatJet_msoftdrop[Matched_FatJet_id]*FatJet_particleNet_massCorr[Matched_FatJet_id]")
data_frame = data_frame.Define("mass_tauBoostedJet_nocorr","FatJet_msoftdrop[Matched_FatJet_id]")


if(int(args.zero_mass) >= 0):
    data_frame = data_frame.Filter("mass_tauBoostedJet_nocorr <= "+args.zero_mass)
    print( "{0} events left in this batch after mass_tauBoostedJet_nocorr ".format( str(data_frame.Count().GetValue()) ) + "<= "+args.zero_mass + " selection.\n")

data_frame = data_frame.Define("Corr_mass","FatJet_particleNet_massCorr[Matched_FatJet_id]")
data_frame = data_frame.Define("TauTauJet_particleNet_XttVsQCD", "FatJet_particleNet_XttVsQCD[Matched_FatJet_id]")


data_frame = data_frame.Define("SubJet_vector", "find_2tau(FatJet_eta, FatJet_phi, SubJet_eta, SubJet_phi, Matched_FatJet_id, SubJet_pt)")
data_frame = data_frame.Define("N_SubJet_Matched", "SubJet_vector[0]")

data_frame = data_frame.Define("boostedtau_vector", "find_2tau(FatJet_eta, FatJet_phi, boostedTau_eta, boostedTau_phi, Matched_FatJet_id, boostedTau_pt)")
data_frame = data_frame.Define("N_boostedtau_Matched", "boostedtau_vector[0]")

data_frame = data_frame.Define("tau_vector", "find_2tau(FatJet_eta, FatJet_phi, Tau_eta, Tau_phi, Matched_FatJet_id, Tau_pt)")
data_frame = data_frame.Define("N_tau_Matched", "tau_vector[0]")

if(args.reco_level_Ana == "y"):

    data_frame = data_frame.Filter("(N_SubJet_Matched > 1) || (N_boostedtau_Matched > 1)");
    print( "{0} events left in this batch after (N_SubJet_Matched > 1) || (N_boostedtau_Matched > 1) selection.\n".format( str(data_frame.Count().GetValue()) ) )
    data_frame = data_frame.Define("mass_CA","M_tautau_CA_boostedtau_SubJet(boostedTau_eta, boostedTau_phi, boostedTau_mass, boostedTau_pt,\
                        SubJet_eta, SubJet_phi, SubJet_mass, SubJet_pt, boostedtau_vector, SubJet_vector, GenMET_pt,GenMET_phi,mass_tauBoostedJet)")


if(args.reco_level_Ana == "n"):
    
    if(args.usingboostedtau == "s"):
        data_frame = data_frame.Filter("N_SubJet_Matched >= 2")
        print( "{0} events left in this batch after N_SubJet_Matched >= 2 selection.\n".format( str(data_frame.Count().GetValue()) ) )

        data_frame = data_frame.Filter("matching_tau_Gentau(GenVisTau_eta, GenVisTau_phi, SubJet_vector, SubJet_eta, SubJet_phi, 0.4) > 0")
        print( "{0} events left in this batch after both_SubJet_Matched_with_GenVistau selection.\n".format( str(data_frame.Count().GetValue()) ) )


    elif(args.usingboostedtau == "b"):
        data_frame = data_frame.Filter("N_boostedtau_Matched >= 2")
        print( "{0} events left in this batch after N_boostedtau_Matched >= 2 selection.\n".format( str(data_frame.Count().GetValue()) ) )

        data_frame = data_frame.Filter("matching_tau_Gentau(GenVisTau_eta, GenVisTau_phi, boostedtau_vector, boostedTau_eta, boostedTau_phi, 0.4) > 0")
        print( "{0} events left in this batch after both_boostedtau_Matched_with_GenVistau selection.\n".format( str(data_frame.Count().GetValue()) ) )


    elif(args.usingboostedtau == "t"):
        data_frame = data_frame.Filter("N_tau_Matched >= 2")
        print( "{0} events left in this batch after N_tau_Matched >= 2 selection.\n".format( str(data_frame.Count().GetValue()) ) )

        data_frame = data_frame.Filter("matching_tau_Gentau(GenVisTau_eta, GenVisTau_phi, tau_vector, Tau_eta, Tau_phi, 0.4) > 0")
        print( "{0} events left in this batch after both_tau_Matched_with_GenVistau selection.\n".format( str(data_frame.Count().GetValue()) ) )


if(args.usingboostedtau == "sb"):
    data_frame = data_frame.Filter("N_SubJet_Matched >= 2")
    print( "{0} events left in this batch after N_SubJet_Matched >= 2 selection.\n".format( str(data_frame.Count().GetValue()) ) )
    
    data_frame = data_frame.Filter("N_boostedtau_Matched >= 2")
    print( "{0} events left in this batch after N_boostedtau_Matched >= 2 selection.\n".format( str(data_frame.Count().GetValue()) ) )

    data_frame = data_frame.Filter("matching_bt_sj(boostedtau_vector, boostedTau_eta, boostedTau_phi, SubJet_vector, SubJet_eta, SubJet_phi, 0.4) > 0")
    print( "{0} events left in this batch after both_boostedtau_Matched_with_subjet selection.\n".format( str(data_frame.Count().GetValue()) ) )


if(args.reco_level_Ana == "n"):
    
    skim_branches =["nGenVisTau","nFatJet","N_Matched_FatJet","N_Matched_1Vistau_FetJet","PNet_ttvsQCD_Matched",

                    "N_SubJet_Matched","N_boostedtau_Matched","N_tau_Matched",
                    
                    "mass_tauBoostedJet","mass_tauBoostedJet_nocorr","Corr_mass"]

    output_file = args.output + "/" + channel + "_" + events_proceed +  "_zeromass_" + args.zero_mass + "_matchtauway_" + args.usingboostedtau + "_1FatJetMatched2VisTau_"+ args.gen_tau_match +"_effiencycheck.root"

elif(args.reco_level_Ana == "y"):
    
    skim_branches =["nGenVisTau","nFatJet", "N_GenTau_Matched",

                    "N_SubJet_Matched","N_boostedtau_Matched","N_tau_Matched",
                    
                    "mass_tauBoostedJet","mass_tauBoostedJet_nocorr","Corr_mass", "mass_CA"]

    output_file = args.output + "/" + channel + "_" + events_proceed + "_recolevel_Ana_" + args.reco_level_Ana + "_zeromass_" + args.zero_mass + "_emuoncut_" + args.emuoncut + "_effiencycheck.root"

data_frame.Snapshot("Events", output_file, skim_branches)

print( "Done")
