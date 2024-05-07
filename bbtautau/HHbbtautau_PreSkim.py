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

ROOT.gInterpreter.Declare('#include "./interface/HHbbtautau_PreSkim.h"')
ROOT.EnableImplicitMT(10)

os.makedirs(args.output, exist_ok=True)
data_base = dasgoclient(args.input)
data_frame = remote_reading([data_base[0]]) 
# data_frame = remote_reading(data_base) 
# Licheng: the problem left here is to separate database into batches when facing a large dataset

print( "We proceed {0} events in this batch\n".format( str(data_frame.Count().GetValue()) ) )

data_frame = data_frame.Define("nGenParticle", "nGenPart")
# number of gen particles

print("Selecting Gen Taus...\n")
data_frame = data_frame.Define("nGenTau", "GenTauHad_Counter(GenPart_pdgId, GenPart_statusFlags)").Filter("nGenTau > 1").Filter("nGenVisTau > 1")
# number of gen taus(events with less than two taus or 2 visibles taus are filtered out)
print( "{0} events left in this batch after Gen and Visible Tau selection.\n".format( str(data_frame.Count().GetValue()) ) )

data_frame = data_frame.Define("SelGenTau_idx", "GenTauHad_Selector(GenPart_pdgId, GenPart_statusFlags)")
# indices of selected taus (this is a vector with length nGenTau)

data_frame = data_frame.Define("SelGenTau_Mother_idx", "GenTauHad_Mother_Selector(GenPart_pdgId, SelGenTau_idx, GenPart_genPartIdxMother, GenPart_statusFlags)")


data_frame = data_frame.Define("SelGenTau_pt_def", "GenTauHad_pt_Def(GenPart_pt, SelGenTau_idx)")
# pt of selected taus in nanoAOD

data_frame = data_frame.Define("GenTauNutrino_idx", "GenTauNutrino_Selector(GenPart_pdgId)")
# indices of selected tau nutrinos (vector)

print("Selecting Gen Tau Nutrinos...\n")
data_frame = data_frame.Define("SelGenTauNutrino_motherIdx", \
                               "SelGenTauNutrino_Mother(GenTauNutrino_idx, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags)")
# indices of the mother of selected tau nutrinos (only if the mother is the tau we need, or the index is -99)
data_frame = data_frame.Define("SelGenTauNutrino_number", "SelGenTauNutrino_Counter(SelGenTauNutrino_motherIdx)").Filter("SelGenTauNutrino_number > 0")
print( "{0} events left in this batch after Gen Tau Nutrinos selection.\n".format( str(data_frame.Count().GetValue()) ) )

print("Selecting Visible Tau matching Gen Tau...\n")
# data_frame = data_frame.Define("GenVisTau_mother_pdgID", \
#                                "GenVisTau_Mother_pdgID(GenVisTau_genPartIdxMother, GenPart_pdgId, GenPart_genPartIdxMother, SelGenTau_Mother_idx, GenPart_statusFlags)")
# data_frame = data_frame.Define("GenVisTau_mother_status", \
#                                "GenVisTau_Mother_status(GenVisTau_genPartIdxMother, GenPart_pdgId, GenPart_genPartIdxMother, SelGenTau_Mother_idx, GenPart_statusFlags)")
data_frame = data_frame.Define("GenVisTau_motherIdx", \
                               "GenVisTau_Mother(GenVisTau_genPartIdxMother, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags)")


data_frame = data_frame.Define("GenVisTau_GenTau_Match", \
                               "GenVisTau_GenTau_Match(GenVisTau_motherIdx, SelGenTau_Mother_idx)").Filter("GenVisTau_GenTau_Match > 1")
data_frame = data_frame.Define("GenTauNutrino_GenTau_Match", \
                               "GenTauNutrino_GenTau_Match(SelGenTauNutrino_motherIdx, SelGenTau_Mother_idx)").Filter("GenTauNutrino_GenTau_Match > 1")
data_frame = data_frame.Define("GenTauNutrino_GenVisTau_Match", \
                               "GenTauNutrino_GenVisTau_Match(SelGenTauNutrino_motherIdx, GenVisTau_motherIdx)").Filter("GenTauNutrino_GenTau_Match > 1")
print( "{0} events left in this batch after Visible Tau/Tau Nutrinos matching Gen Tau selection.\n".format( str(data_frame.Count().GetValue()) ) )

data_frame = data_frame.Define("Leading_pT_GenVisTau_Idx", "Leading_pT_GenVisTau_Idx_selecter(GenVisTau_pt)")

data_frame = data_frame.Define("Mass_HTauTau_GenVisTau", "Mass_HTauTau_GenVisTau_Calc(GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, Leading_pT_GenVisTau_Idx)")

data_frame = data_frame.Define("Mass_HTauTau_CA", "Mass_HTauTau_CA_calc(Mass_HTauTau_GenVisTau, Leading_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenVisTau_pt, GenPart_pt)")


print( "{0} events left in this batch.\n".format( str(data_frame.Count().GetValue()) ) )

skim_branches = ["nGenParticle", "nGenTau", "nGenVisTau", "SelGenTauNutrino_number",
                 
                 "SelGenTau_idx", "GenTauNutrino_idx", "SelGenTau_Mother_idx", "SelGenTauNutrino_motherIdx", "GenVisTau_motherIdx", 
                 
                 "GenVisTau_pt", "GenVisTau_eta", "GenVisTau_phi", "GenVisTau_mass",
                 "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass",
                 
                 "GenJet_mass", "GenJet_pt", "GenJet_eta", "GenJet_phi", "GenJet_partonFlavour",
                 
                 "SelGenTau_pt_def",
                 
                 "Leading_pT_GenVisTau_Idx",
                 
                 "Mass_HTauTau_GenVisTau", "Mass_HTauTau_CA"
                ]

output_file = args.output + "/" + "HHbbtautau_PreSkim.root"
data_frame.Snapshot("Events", output_file, skim_branches)

