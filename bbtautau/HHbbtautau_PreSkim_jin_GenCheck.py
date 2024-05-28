import os, sys
import subprocess
import ROOT
import argparse
import time
#from utils.CMShelper import *

parser = argparse.ArgumentParser(description='Skim full tuple.')
parser.add_argument('-i', '--input', required=False, type=str, default="/VBFHHto2B2Tau_CV-1_C2V-1_C3-1_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM", help="output directory name")
parser.add_argument('-o', '--output', required=False, type=str, default="./OutPut", help="output directory name")
args = parser.parse_args()

ROOT.gInterpreter.Declare('#include "./interface/HHbbtautau_PreSkim_jin_Gen.h"')
ROOT.EnableImplicitMT(10)

os.makedirs(args.output, exist_ok=True)
#data_base = dasgoclient(args.input)
#data_frame = remote_reading([data_base[0]]) 


data_frame = ROOT.RDataFrame('Events', "/afs/cern.ch/user/l/lichengz/public/forbbtautau/294beabb-42f5-4d00-898c-9ed7e0e90049.root")
###this file error
#data_frame = ROOT.RDataFrame('Events', "/afs/cern.ch/user/l/lichengz/public/forbbtautau/48b43825-6a4c-49f1-836d-f382204dc607.root")
#data_frame = ROOT.RDataFrame('Events', "/afs/cern.ch/user/l/lichengz/public/forbbtautau/4fea44a0-387e-453d-b87c-3619a9c6ae5d.root")
#data_frame = ROOT.RDataFrame('Events', "/afs/cern.ch/user/l/lichengz/public/forbbtautau/b520ac03-8ab1-4321-8819-581106bef67c.root")


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

data_frame = data_frame.Define("GenTauNutrino_idx", "GenTauNutrino_Selector(GenPart_pdgId, GenPart_statusFlags)")
# indices of selected tau nutrinos (vector)

print("Selecting Gen Tau Nutrinos...\n")
data_frame = data_frame.Define("SelGenTauNutrino_motherIdx", \
                               "SelGenTauNutrino_Mother(GenTauNutrino_idx, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags)")

data_frame = data_frame.Define("SelGenTauNutrino_pT", \
                               "SelGenTauNutrino_pT(SelGenTauNutrino_motherIdx, GenPart_pt)")
# indices of the mother of selected tau nutrinos (only if the mother is the tau we need, or the index is -99)
data_frame = data_frame.Define("SelGenTauNutrino_number", "SelGenTauNutrino_Counter(SelGenTauNutrino_motherIdx)").Filter("SelGenTauNutrino_number > 0")
print( "{0} events left in this batch after Gen Tau Nutrinos selection.\n".format( str(data_frame.Count().GetValue()) ) )

print("Selecting Visible Tau matching Gen Tau...\n")
## data_frame = data_frame.Define("GenVisTau_mother_pdgID", \
#                                "GenVisTau_Mother_pdgID(GenVisTau_genPartIdxMother, GenPart_pdgId, GenPart_genPartIdxMother, SelGenTau_Mother_idx, GenPart_statusFlags)")
## data_frame = data_frame.Define("GenVisTau_mother_status", \
#                                "GenVisTau_Mother_status(GenVisTau_genPartIdxMother, GenPart_pdgId, GenPart_genPartIdxMother, SelGenTau_Mother_idx, GenPart_statusFlags)")
data_frame = data_frame.Define("GenVisTau_motherIdx", \
                               "GenVisTau_Mother(GenVisTau_genPartIdxMother, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags)")


data_frame = data_frame.Define("GenVisTau_GenTau_Match", \
                               "GenVisTau_GenTau_Match(GenVisTau_motherIdx, SelGenTau_Mother_idx)").Filter("GenVisTau_GenTau_Match > 1")
data_frame = data_frame.Define("GenTauNutrino_GenTau_Match", \
                               "GenTauNutrino_GenTau_Match(SelGenTauNutrino_motherIdx, SelGenTau_Mother_idx)").Filter("GenTauNutrino_GenTau_Match > 1")
data_frame = data_frame.Define("GenTauNutrino_GenVisTau_Match", "GenTauNutrino_GenVisTau_Match(SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_motherIdx, GenPart_pt)").Filter("GenTauNutrino_GenVisTau_Match > 1")
print( "{0} events left in this batch after Visible Tau/Tau Nutrinos matching Gen Tau selection.\n".format( str(data_frame.Count().GetValue()) ) )

###Leading_pt_visTau
####data_frame = data_frame.Define("Leading_pT_GenVisTau_Idx", "Leading_pT_GenVisTau_Idx_selecter(GenVisTau_pt)")

####data_frame = data_frame.Define("Mass_HTauTau_GenVisTau", "Mass_HTauTau_GenVisTau_Calc(GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, Leading_pT_GenVisTau_Idx)")

####data_frame = data_frame.Define("Mass_HTauTau_CA", "Mass_HTauTau_CA_calc(Mass_HTauTau_GenVisTau, Leading_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenPart_pt)")

####data_frame = data_frame.Define("Mass_HTauTau_CA_all_N", "Mass_HTauTau_CA_calc_all_N(Mass_HTauTau_GenVisTau, Leading_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")


data_frame = data_frame.Define("SameHMother_GenVisTau_Idx", "SameHMother_GenVisTau_Idx_selecter(GenVisTau_motherIdx, GenPart_pdgId, GenPart_genPartIdxMother)")

data_frame = data_frame.Define("SameHMother_pT_GenVisTau_Idx", "SameHMother_pT_GenVisTau_Idx_selecter(GenVisTau_motherIdx, GenPart_pdgId, GenPart_genPartIdxMother)").Filter("SameHMother_GenVisTau_Idx > 1")

data_frame = data_frame.Define("Mass_HTauTau_GenVisTau_H_mother", "Mass_HTauTau_GenVisTau_Calc(GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, SameHMother_pT_GenVisTau_Idx)")

##for each tau, only one Nu used
####data_frame = data_frame.Define("Mass_HTauTau_CA_H_mother", "Mass_HTauTau_CA_calc(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenPart_pt)")


##MH_CA
##data_frame = data_frame.Define("Mass_HTauTau_CA_all_N_H_mother", "Mass_HTauTau_CA_calc_all_N(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
data_frame = data_frame.Define("Mass_HTauTau_CA_all_N_H_mother_CACA", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[19]")
data_frame = data_frame.Define("Mass_HTauTau_CA_METtoNu", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[20]")
data_frame = data_frame.Define("Mass_HTauTau_CA_METtoNu_tau", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[21]")



print( "{0} events left in this batch after Gen and Visible Tau pt selection.\n".format( str(data_frame.Count().GetValue()) ) )

##MH_4V
data_frame = data_frame.Define("Mass_HTauTau_4V_all_N_H_mother_plus4V", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[0]")

##dr & pt check
data_frame = data_frame.Define("dr_vistau_nutrino_0", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[1]")

data_frame = data_frame.Define("dr_vistau_nutrino_1", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[2]")

data_frame = data_frame.Define("pt_vistau_0", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[3]")

data_frame = data_frame.Define("pt_vistau_1", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[4]")

data_frame = data_frame.Define("pt_H_CA", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[5]")

data_frame = data_frame.Define("dr_vistau01", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[6]")


#Nu check
data_frame = data_frame.Define("Pt_Nu0", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[7]")

data_frame = data_frame.Define("Pt_Nu1", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[8]")

data_frame = data_frame.Define("M_Nu0", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[9]")

data_frame = data_frame.Define("M_Nu1", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[10]")


#MET to Nu
data_frame = data_frame.Define("Pt_MET_Nu0", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[11]")
data_frame = data_frame.Define("d_Pt_MET_Nu0", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[15]")
data_frame = data_frame.Define("Pt_MET_Nu1", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[12]")
data_frame = data_frame.Define("d_Pt_MET_Nu1", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[16]")
data_frame = data_frame.Define("Pt_MET_Nu0_tau", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[13]")
data_frame = data_frame.Define("d_Pt_MET_Nu0_tau", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[17]")
data_frame = data_frame.Define("Pt_MET_Nu1_tau", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[14]")
data_frame = data_frame.Define("d_Pt_MET_Nu1_tau", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[18]")


data_frame = data_frame.Define("d_Phi_MET_Nu0", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[22]")

data_frame = data_frame.Define("d_Phi_MET_Nu1", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[23]")

data_frame = data_frame.Define("d_Phi_Nu1_Nu0", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[24]")



data_frame = data_frame.Define("MET_from_Nu01_Pt", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[25]")

data_frame = data_frame.Define("MET_from_Nu01_Phi", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[26]")


data_frame = data_frame.Define("d_Pt_GenMET_NuMET", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[27]")

data_frame = data_frame.Define("d_Phi_GenMET_NuMET", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[28]")


##all_Nu
data_frame = data_frame.Define("Pt_All_Nu", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[29]")
data_frame = data_frame.Define("Phi_All_Nu", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[30]")

data_frame = data_frame.Define("Number_NutoTau0", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[31]")
data_frame = data_frame.Define("Number_NutoTau1", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[32]")
data_frame = data_frame.Define("Number_AllNu", \
                               "Mass_HTauTau_CA_calc_all_N_plus4V(Mass_HTauTau_GenVisTau_H_mother, SameHMother_pT_GenVisTau_Idx, GenVisTau_motherIdx, SelGenTauNutrino_motherIdx, GenTauNutrino_idx, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenMET_pt, GenMET_phi)[33]")

#pt_cut
#data_frame = data_frame.Define("GenVisTau_pT_cut", "pT_min_cut(GenVisTau_pt)").Filter("pt_vistau_1 > (float)30.0").Filter("pt_vistau_0 > (float)30.0")
#print( "{0} events left in this batch after Gen and Visible Tau pt selection.\n".format( str(data_frame.Count().GetValue()) ) )


print( "{0} events left in this batch.\n".format( str(data_frame.Count().GetValue()) ) )

'''skim_branches = ["nGenParticle", "nGenTau", "nGenVisTau", "SelGenTauNutrino_number",
                 
                 "SelGenTau_idx", "GenTauNutrino_idx", "SelGenTau_Mother_idx", "SelGenTauNutrino_motherIdx", "GenVisTau_motherIdx", 
                 
                 "GenVisTau_pt", "GenVisTau_eta", "GenVisTau_phi", "GenVisTau_mass",
                 "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass",
                 
                 "GenJet_mass", "GenJet_pt", "GenJet_eta", "GenJet_phi", "GenJet_partonFlavour",
                 
                 "SelGenTau_pt_def",
                 
                 "Leading_pT_GenVisTau_Idx",
                 
                 "Mass_HTauTau_GenVisTau", "Mass_HTauTau_CA"
                ]'''
skim_branches = ["Mass_HTauTau_GenVisTau_H_mother", "Mass_HTauTau_CA_all_N_H_mother_CACA", "Mass_HTauTau_4V_all_N_H_mother_plus4V", 

                "dr_vistau_nutrino_0", "dr_vistau_nutrino_1","pt_vistau_0","pt_vistau_1","pt_H_CA","dr_vistau01",
                
                "Pt_Nu0","Pt_Nu1","M_Nu0","M_Nu1",
                
                "GenMET_phi","GenMET_pt", "Pt_MET_Nu0", "d_Pt_MET_Nu0", "Pt_MET_Nu1", "d_Pt_MET_Nu1",
                
                "Pt_MET_Nu0_tau", "d_Pt_MET_Nu0_tau", "Pt_MET_Nu1_tau", "d_Pt_MET_Nu1_tau",
                
                "Mass_HTauTau_CA_METtoNu","Mass_HTauTau_CA_METtoNu_tau",
                
                "d_Phi_MET_Nu0", "d_Phi_MET_Nu1", "d_Phi_Nu1_Nu0", "MET_from_Nu01_Pt", "MET_from_Nu01_Phi", "d_Pt_GenMET_NuMET", "d_Phi_GenMET_NuMET",
                
                "Number_NutoTau0", "Number_NutoTau1", "Number_AllNu", "Pt_All_Nu", "Phi_All_Nu",
                
                "boostedTau_mass", "boostedTau_phi", "boostedTau_eta", "boostedTau_pt",
                
                "FatJet_eta", "FatJet_mass", "FatJet_phi", "FatJet_pt"]

output_file = args.output + "/" + "HHbbtautau_PreSkim_Gen.root"
data_frame.Snapshot("Events", output_file, skim_branches)

