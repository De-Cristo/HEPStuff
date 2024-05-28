#include "ROOT/RVec.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"
#include <nlohmann/json.hpp>
#include <iostream>
#include <vector>
#include <algorithm> // for std::max_element
#include <cmath>

using namespace ROOT;
using namespace ROOT::VecOps;
using RNode = ROOT::RDF::RNode;
using str = const std::string &;

using cRVecF = const ROOT::RVecF &;
using cRVecI = const ROOT::RVecI &;
using cRVecC = const ROOT::RVecC &;
using cRVecU = const ROOT::RVecU &;
using cRVecB = const ROOT::RVecB &;

// std::vector<int> myVector;

bool check_bit(int number, int bitpos) {
    int res = number & (1 << bitpos);
    return res != 0;
}


bool isNumberInList(RVecI myList, int target) {
    for (int i = 0; i < myList.size(); ++i) {
        if (myList[i] == target) {
            return true;
        }
    }
    return false;
}

std::pair<int, int> findLeadingAndSubleading(cRVecF vec) {
    // Initialize variables to store the indices of leading and subleading elements
    int leadingIdx = 0;
    int subleadingIdx = 0;

    // Find the index of the leading element
    auto leadingIter = std::max_element(vec.begin(), vec.end());
    leadingIdx = std::distance(vec.begin(), leadingIter);

    // Find the index of the subleading element
    int leadingValue = *leadingIter;
    auto subleadingIter = std::max_element(vec.begin(), vec.end(), [&](int a, int b) {
        return a < b && b < leadingValue;
    });
    subleadingIdx = std::distance(vec.begin(), subleadingIter);

    // Return a pair containing the indices of leading and subleading elements
    return std::make_pair(leadingIdx, subleadingIdx);
}


std::vector<int> findIndicesOfValue(cRVecF vec, int value) {
    std::vector<int> indices;
    for (int i = 0; i < vec.size(); ++i) {
        if (vec[i] == value) {
            indices.push_back(i);
        }
    }
    return indices;
}

float deltaR2(float eta_1, float eta_2, float phi_1, float phi_2){
    const float deta = eta_1 - eta_2;
    const float dphi = ROOT::Math::VectorUtil::Phi_mpi_pi(phi_1 - phi_2);
    const float dRsq = std::pow(deta,2) + std::pow(dphi,2);
    return dRsq;
}

float deltaR(float eta_1, float eta_2, float phi_1, float phi_2){
    const float deta = eta_1 - eta_2;
    const float dphi = ROOT::Math::VectorUtil::Phi_mpi_pi(phi_1 - phi_2);
    const float dRsq = std::pow(deta,2) + std::pow(dphi,2);
    return sqrt(dRsq);
}


int check_tau_mother(cRVecI GenPart_pdgId, cRVecI GenPart_genPartIdxMother, cRVecU GenPart_statusFlags, int idx, string isTau) {
    if (isTau == "Tau") {
        if ( (check_bit(GenPart_statusFlags[idx], 12)) ) {
            return idx;
        }
        if ( (abs(GenPart_pdgId[GenPart_genPartIdxMother[idx]]) == 15) && (check_bit(GenPart_statusFlags[GenPart_genPartIdxMother[idx]], 12)) ) { 
            return GenPart_genPartIdxMother[idx]; 
        }
        
        return check_tau_mother(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags, GenPart_genPartIdxMother[idx], isTau);
    }
    
    if (isTau == "NuTau") {
        if ( (abs(GenPart_pdgId[GenPart_genPartIdxMother[idx]]) == 15) && (check_bit(GenPart_statusFlags[GenPart_genPartIdxMother[idx]], 12)) ) { 
            return GenPart_genPartIdxMother[idx];
        }
        if ( (abs(GenPart_pdgId[GenPart_genPartIdxMother[idx]]) != 15) ) {
            return -99;
        }
        
        return check_tau_mother(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags, GenPart_genPartIdxMother[idx], isTau);
    }
    
    if (isTau == "visTau") {
        if ( (abs(GenPart_pdgId[idx]) == 15) && (check_bit(GenPart_statusFlags[idx], 12)) ) {
            return idx;
        }
        if ( (abs(GenPart_pdgId[GenPart_genPartIdxMother[idx]]) == 15) && (check_bit(GenPart_statusFlags[GenPart_genPartIdxMother[idx]], 12)) ) {
            return GenPart_genPartIdxMother[idx];
        }
        if ( (abs(GenPart_pdgId[GenPart_genPartIdxMother[idx]]) != 15) ) {
            return idx;
        }
        
        return check_tau_mother(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags, GenPart_genPartIdxMother[idx], isTau);
    }
    
    return -99;
}


int GenTauHad_Counter(cRVecI GenPart_pdgId,
                      cRVecU GenPart_statusFlags){
    
    int _nGentau = 0;
        
    for(auto i = 0; i < GenPart_pdgId.size(); i++) {
        
        if ( (abs(GenPart_pdgId[i]) == 15) && (check_bit(GenPart_statusFlags[i], 13)) ) {
            _nGentau++;
        }
    }

    return _nGentau;
}


RVecI GenTauHad_Selector( cRVecI GenPart_pdgId, 
                          cRVecU GenPart_statusFlags){
    
    RVecI GenTauIdx;
    int idx = -99;
    for(auto i = 0; i < GenPart_pdgId.size(); i++) {
        
        if ( (abs(GenPart_pdgId[i]) == 15) && (check_bit(GenPart_statusFlags[i], 13)) ) {
             idx = i;
             GenTauIdx.push_back(idx);
        }
        
    }

    return GenTauIdx; // GenTauIdx(idx_tau1, idx_tau2, ...)
}


RVecI GenTauHad_Mother_Selector(cRVecI GenPart_pdgId,
                                cRVecI SelGenTau_idx, 
                                cRVecI GenPart_genPartIdxMother, 
                                cRVecU GenPart_statusFlags){

    RVecI GenTau_motherIdx;
    int idx = -99;
    str isTau = "Tau";
    for(auto i = 0; i < SelGenTau_idx.size(); i++) {
        idx = check_tau_mother(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags, SelGenTau_idx[i], isTau);
        GenTau_motherIdx.push_back(idx);
    }

    return GenTau_motherIdx;
}


RVecF GenTauHad_pt_Def(RVecF GenPart_pt, 
                       RVecI SelGenTau_idx){
    RVecF GenTauPt;
    
    for (auto i = 0; i < SelGenTau_idx.size(); i++) {
        GenTauPt.push_back(GenPart_pt[SelGenTau_idx[i]]);
    }
    
    return GenTauPt;
}


RVecI GenTauNutrino_Selector(cRVecI GenPart_pdgId,
                             cRVecU GenPart_statusFlags){

    RVecI GenNutrinoIdx;
    
    int idx = -99;
    for(auto i = 0; i < GenPart_pdgId.size(); i++) {
        if ( (abs(GenPart_pdgId[i]) == 16) && (check_bit(GenPart_statusFlags[i], 13)) ) {
            idx = i;
            GenNutrinoIdx.push_back(idx);
        }
    }

    return GenNutrinoIdx; // GenTauIdx(idx_nvtau1, idx_nvtau2, ...)
}


RVecI SelGenTauNutrino_Mother(cRVecI GenTauNutrino_idx , 
                              cRVecI GenPart_pdgId,
                              cRVecI GenPart_genPartIdxMother, 
                              cRVecU GenPart_statusFlags ){
    RVecI GenNutrino_motherIdx;
    int idx = -99;
    str isTau = "NuTau";
    for(auto i = 0; i < GenTauNutrino_idx.size(); i++) {
        idx = check_tau_mother(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags, GenTauNutrino_idx[i], isTau);
        GenNutrino_motherIdx.push_back(idx);
    }
    
    return GenNutrino_motherIdx;
}

RVecF SelGenTauNutrino_pT(cRVecI SelGenTauNutrino_motherIdx , cRVecF GenPart_pt ){
    RVecF GenNutrino_pT;
    for(auto i = 0; i < SelGenTauNutrino_motherIdx.size(); i++) {
        if (SelGenTauNutrino_motherIdx[i] != -99){
            GenNutrino_pT.push_back(GenPart_pt[i]);
            // std::cout<< "with TauMother: "<< GenPart_pt[i]<<std::endl;
        }
        else{
            GenNutrino_pT.push_back(-GenPart_pt[i]);
            // std::cout<<GenPart_pt[i]<<std::endl;
        }
    }
    return GenNutrino_pT;
}


// RVecI GenVisTau_Mother_pdgID(cRVecF GenVisTau_genPartIdxMother,
//                        cRVecI GenPart_pdgId,
//                        cRVecF GenPart_genPartIdxMother, 
//                        cRVecI SelGenTau_Mother_idx,
//                        cRVecU GenPart_statusFlags){

//     RVecI GenVisTau_mother_pdgID;
//     for (auto i = 0; i < GenVisTau_genPartIdxMother.size(); i++) {
//         GenVisTau_mother_pdgID.push_back(GenPart_pdgId[GenVisTau_genPartIdxMother[i]]);
//     }
    
//     return GenVisTau_mother_pdgID;
// }


// RVecI GenVisTau_Mother_status(cRVecF GenVisTau_genPartIdxMother,
//                        cRVecI GenPart_pdgId,
//                        cRVecF GenPart_genPartIdxMother, 
//                        cRVecI SelGenTau_Mother_idx,
//                        cRVecU GenPart_statusFlags){

//     RVecI GenVisTau_mother_status;
    
//     for (auto i = 0; i < GenVisTau_genPartIdxMother.size(); i++) {
//         if ( check_bit(GenPart_statusFlags[GenVisTau_genPartIdxMother[i]], 12) ) {
//             GenVisTau_mother_status.push_back(0);
//         }else{
//             if ( check_bit(GenPart_statusFlags[GenVisTau_genPartIdxMother[i]], 13) ) {
//                 GenVisTau_mother_status.push_back(13);
//             }else{
//                 if ( check_bit(GenPart_statusFlags[GenVisTau_genPartIdxMother[i]], 7) ) {
//                     GenVisTau_mother_status.push_back(7);
//                 }else{
//                     GenVisTau_mother_status.push_back(15);
//                 }
//             }
//         }
//     }
    
//     return GenVisTau_mother_status;
// }


RVecI GenVisTau_Mother(cRVecI GenVisTau_genPartIdxMother,
                       cRVecI GenPart_pdgId,
                       cRVecI GenPart_genPartIdxMother, 
                       cRVecU GenPart_statusFlags){

    RVecI GenVisTau_motherIdx;
    int idx = -99;
    str isTau = "visTau";
    for(auto i = 0; i < GenVisTau_genPartIdxMother.size(); i++) {
        idx = check_tau_mother(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags, GenVisTau_genPartIdxMother[i], isTau);
        GenVisTau_motherIdx.push_back(idx);
    }
    
    return GenVisTau_motherIdx;
}


int SelGenTauNutrino_Counter(cRVecI SelGenNutrino_motherIdx){
    int _nSelTauNutrino = 0;
    
    //std::cout << "pT of matched Nu = " << std::endl;
    
    for(auto i = 0; i < SelGenNutrino_motherIdx.size(); i++) {
        if ( SelGenNutrino_motherIdx[i] != -99 ) { _nSelTauNutrino++; }
    }
    
    return _nSelTauNutrino;
}

int GenVisTau_GenTau_Match(cRVecI GenVisTau_motherIdx, 
                           cRVecI SelGenTau_Mother_idx){

    int _nMatch = 0;

    for (auto i = 0; i < GenVisTau_motherIdx.size(); i++) {
        if (isNumberInList(SelGenTau_Mother_idx, GenVisTau_motherIdx[i])) {
            _nMatch++;
        }
    }
    
    return _nMatch;
        
}

int GenTauNutrino_GenTau_Match(cRVecI SelGenTauNutrino_motherIdx, 
                               cRVecI SelGenTau_Mother_idx){

    //std::cout << "pT of matched Nu = " << std::endl;
    
    int _nMatch = 0;


    for (auto i = 0; i < SelGenTauNutrino_motherIdx.size(); i++) {
        if (isNumberInList(SelGenTau_Mother_idx, SelGenTauNutrino_motherIdx[i])) {
            _nMatch++;
        }
    }
    
    return _nMatch;
        
}

int GenTauNutrino_GenVisTau_Match(cRVecI SelGenTauNutrino_motherIdx, 
                                  cRVecI GenTauNutrino_idx,
                                  cRVecI GenVisTau_motherIdx,
                                  cRVecF GenPart_pt){

    int _nMatch = 0;
    for (auto i = 0; i < SelGenTauNutrino_motherIdx.size(); i++) {
        if (isNumberInList(GenVisTau_motherIdx, SelGenTauNutrino_motherIdx[i])) {
            _nMatch++;
            //std::cout<<"pT of matched nu = " << GenPart_pt[GenTauNutrino_idx[i]] << std::endl;
        }
    }
    
    return _nMatch;
}

RVecI Leading_pT_GenVisTau_Idx_selecter(cRVecF GenVisTau_pt){

    RVecI GenVisTau_Idx;
    
    std::pair<int, int> temp_idx_pair;
    temp_idx_pair = findLeadingAndSubleading(GenVisTau_pt);
    
    GenVisTau_Idx.push_back(temp_idx_pair.first);
    GenVisTau_Idx.push_back(temp_idx_pair.second);
    
    return GenVisTau_Idx;
}

int check_Higgs_mother(int idx,
                       cRVecI GenPart_pdgId,
                       cRVecI GenPart_genPartIdxMother) {

    if ( abs(GenPart_pdgId[GenPart_genPartIdxMother[idx]]) == 25 ) {
        return GenPart_genPartIdxMother[idx];
    }
    else {
        return -abs(GenPart_pdgId[GenPart_genPartIdxMother[idx]]);
    }
}


int SameHMother_GenVisTau_Idx_selecter(cRVecI GenVisTau_motherIdx, cRVecI GenPart_pdgId,
                       cRVecI GenPart_genPartIdxMother){

    int H_mother_number = 0;
    std::vector<int> H_idx_vec ;
    
    for(int idx : GenVisTau_motherIdx) {
        int H_idx;
        H_idx = check_Higgs_mother(idx, GenPart_pdgId, GenPart_genPartIdxMother);
        if(H_idx >= 0) {
            H_mother_number += 1;
            H_idx_vec.push_back(H_idx);
        }
    }
    if (H_idx_vec.size() == 2){
        if (H_idx_vec[0] != H_idx_vec[1]) {
            H_mother_number = -H_mother_number;
        }
    }

    return H_mother_number;
}

RVecI SameHMother_pT_GenVisTau_Idx_selecter(cRVecI GenVisTau_motherIdx, cRVecI GenPart_pdgId,
                       cRVecI GenPart_genPartIdxMother){

    std::vector<int> Tau_idx_vec ;
    
    for (int i = 0; i < GenVisTau_motherIdx.size() ; ++i) {
        int H_idx;
        H_idx = check_Higgs_mother(GenVisTau_motherIdx[i], GenPart_pdgId, GenPart_genPartIdxMother);
        if(H_idx >= 0) {
            Tau_idx_vec.push_back(i);
        }
    }
    

    return Tau_idx_vec ;
}


float Mass_HTauTau_GenVisTau_Calc(cRVecF GenVisTau_pt,
                                  cRVecF GenVisTau_eta,
                                  cRVecF GenVisTau_phi,
                                  cRVecF GenVisTau_mass,
                                  cRVecI Leading_pT_GenVisTau_Idx){
    float _Mass_HTauTau_GenVisTau = 0;
    
    ROOT::Math::PtEtaPhiMVector v4_vis1(GenVisTau_pt[Leading_pT_GenVisTau_Idx[0]], 
                       GenVisTau_eta[Leading_pT_GenVisTau_Idx[0]], 
                       GenVisTau_phi[Leading_pT_GenVisTau_Idx[0]], 
                       GenVisTau_mass[Leading_pT_GenVisTau_Idx[0]]);
    
    ROOT::Math::PtEtaPhiMVector v4_vis2(GenVisTau_pt[Leading_pT_GenVisTau_Idx[1]],
                       GenVisTau_eta[Leading_pT_GenVisTau_Idx[1]],
                       GenVisTau_phi[Leading_pT_GenVisTau_Idx[1]],
                       GenVisTau_mass[Leading_pT_GenVisTau_Idx[1]]);
    
    ROOT::Math::PtEtaPhiMVector v4_HTauTau = v4_vis1 + v4_vis2;
    
    _Mass_HTauTau_GenVisTau = v4_HTauTau.M();
    
    return _Mass_HTauTau_GenVisTau;
}



void DecomposeMomentum(
    const ROOT::Math::PtEtaPhiMVector& p, 
    const ROOT::Math::PtEtaPhiMVector& a, 
    const ROOT::Math::PtEtaPhiMVector& b,
    ROOT::Math::PtEtaPhiMVector& pa,
    ROOT::Math::PtEtaPhiMVector& pb
) {
    //
    typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> CartesianVector;

    //
    CartesianVector p_xyz(p.Px(), p.Py(), p.Pz(), p.E());
    CartesianVector a_xyz(a.Px(), a.Py(), a.Pz(), a.E());
    CartesianVector b_xyz(b.Px(), b.Py(), b.Pz(), b.E());

    //
    double a_magnitude = a_xyz.R();
    double b_magnitude = b_xyz.R();

    // 
    CartesianVector a_unit = a_xyz / a_magnitude;
    CartesianVector b_unit = b_xyz / b_magnitude;

    //
    double p_dot_a = p_xyz.Dot(a_unit);
    double p_dot_b = p_xyz.Dot(b_unit);

    CartesianVector p_a = p_dot_a * a_unit;
    CartesianVector p_b = p_dot_b * b_unit;

    //
    pa = ROOT::Math::PtEtaPhiMVector(p_a.Pt(), p_a.Eta(), p_a.Phi(), p_a.M());
    pb = ROOT::Math::PtEtaPhiMVector(p_b.Pt(), p_b.Eta(), p_b.Phi(), p_b.M());
}


std::vector<float> Mass_HTauTau_CA_calc_all_N_plus4V(float Mass_HTauTau_GenVisTau,
                           cRVecI Leading_pT_GenVisTau_Idx,
                           cRVecI GenVisTau_motherIdx,
                           cRVecI SelGenTauNutrino_motherIdx,
                           cRVecI GenNutrinoIdx,
                           cRVecF GenVisTau_pt,
                           cRVecF GenVisTau_eta,
                           cRVecF GenVisTau_phi,
                           cRVecF GenVisTau_mass,
                           cRVecF GenPart_pt,
                           cRVecF GenPart_eta,
                           cRVecF GenPart_phi,
                           cRVecF GenPart_mass,
                           float GenMET_pt,
                           float GenMET_phi
                           //cRVecI GenPart_pdgId,
                           //cRVecU GenPart_statusFlags
                           ){
    
    std::vector<float> mess_dr_01;
    float _Mass_HTauTau_4V = 0;
    float dr_0 = 0;
    float dr_1 = 0;
    
    int mother0_idx = GenVisTau_motherIdx[Leading_pT_GenVisTau_Idx[0]];
    int mother1_idx = GenVisTau_motherIdx[Leading_pT_GenVisTau_Idx[1]];
    
    std::vector<int> NuTau_indices_forM0 = findIndicesOfValue(SelGenTauNutrino_motherIdx, mother0_idx);
    std::vector<int> NuTau_indices_forM1 = findIndicesOfValue(SelGenTauNutrino_motherIdx, mother1_idx);

    //std::cout<< "matched Nutau0 number: " << NuTau_indices_forM0.size() << std::endl;
    //std::cout<< "matched Nutau1 number: " << NuTau_indices_forM1.size() << std::endl;
    
    ROOT::Math::PtEtaPhiMVector v4_AllN_forM0;
    ROOT::Math::PtEtaPhiMVector v4_AllN_forM1;

    for (int i : NuTau_indices_forM0) {
        int idx = GenNutrinoIdx[i];
        ROOT::Math::PtEtaPhiMVector v4_N(GenPart_pt[idx], 
                       GenPart_eta[idx], 
                       GenPart_phi[idx], 
                       GenPart_mass[idx]);
        
        v4_AllN_forM0 += v4_N;
    }


    for (int i : NuTau_indices_forM1) {
        int idx = GenNutrinoIdx[i];
        ROOT::Math::PtEtaPhiMVector v4_N(GenPart_pt[idx], 
                       GenPart_eta[idx], 
                       GenPart_phi[idx], 
                       GenPart_mass[idx]);
        
        v4_AllN_forM1 += v4_N;
    }

    ROOT::Math::PtEtaPhiMVector v4_tau_forM0(GenVisTau_pt[Leading_pT_GenVisTau_Idx[0]], GenVisTau_eta[Leading_pT_GenVisTau_Idx[0]], 
                             GenVisTau_phi[Leading_pT_GenVisTau_Idx[0]], GenVisTau_mass[Leading_pT_GenVisTau_Idx[0]]);

    ROOT::Math::PtEtaPhiMVector v4_tau_forM1(GenVisTau_pt[Leading_pT_GenVisTau_Idx[1]], GenVisTau_eta[Leading_pT_GenVisTau_Idx[1]], 
                             GenVisTau_phi[Leading_pT_GenVisTau_Idx[1]], GenVisTau_mass[Leading_pT_GenVisTau_Idx[1]]);

    ROOT::Math::PtEtaPhiMVector v4_HTauTau_CA = v4_AllN_forM0 + v4_AllN_forM1 + v4_tau_forM0 + v4_tau_forM1;
    
    _Mass_HTauTau_4V = v4_HTauTau_CA.M();

    // dr_0 = std::sqrt((v4_AllN_forM0.Phi()-v4_tau_forM0.Phi())*(v4_AllN_forM0.Phi()-v4_tau_forM0.Phi())-(v4_AllN_forM0.Eta()-v4_tau_forM0.Eta())*(v4_AllN_forM0.Eta()-v4_tau_forM0.Eta()));
    
    dr_0 = deltaR(v4_AllN_forM0.Eta(), v4_tau_forM0.Eta(), v4_AllN_forM0.Phi(), v4_tau_forM0.Phi());
    // float deltaR(float eta_1, float eta_2, float phi_1, float phi_2){
    
    // dr_1 = std::sqrt((v4_AllN_forM1.Phi()-v4_tau_forM1.Phi())*(v4_AllN_forM1.Phi()-v4_tau_forM1.Phi())-(v4_AllN_forM1.Eta()-v4_tau_forM1.Eta())*(v4_AllN_forM1.Eta()-v4_tau_forM1.Eta()));
    dr_1 = deltaR(v4_AllN_forM1.Eta(), v4_tau_forM1.Eta(), v4_AllN_forM1.Phi(), v4_tau_forM1.Phi());

    float _Pt_HTauTau_CA = v4_HTauTau_CA.Pt();
    float dr_01 = deltaR(v4_tau_forM1.Eta(), v4_tau_forM0.Eta(), v4_tau_forM1.Phi(), v4_tau_forM0.Phi());


    mess_dr_01.push_back(_Mass_HTauTau_4V);

    mess_dr_01.push_back(dr_0);
    mess_dr_01.push_back(dr_1);
    mess_dr_01.push_back(v4_tau_forM0.Pt());
    mess_dr_01.push_back(v4_tau_forM1.Pt());

    mess_dr_01.push_back(_Pt_HTauTau_CA);
    mess_dr_01.push_back(dr_01);

    mess_dr_01.push_back(v4_AllN_forM0.Pt());
    mess_dr_01.push_back(v4_AllN_forM1.Pt());

    mess_dr_01.push_back(v4_AllN_forM0.M());
    mess_dr_01.push_back(v4_AllN_forM1.M());


    //MET to Nu

    float pt_MET = GenMET_pt;
    float phi_MET = GenMET_phi;
    float Nu0_pt = 0;
    float Nu1_pt = 0;

    float dphi0;
    float dphi1;
    float dphi01;

    ROOT::Math::PtEtaPhiMVector p(pt_MET, 0, phi_MET, pt_MET);
    ROOT::Math::PtEtaPhiMVector pa, pb;

    DecomposeMomentum(p, v4_AllN_forM0, v4_AllN_forM1, pa, pb);

    //Nu0_pt = std::fabs(pt_MET / std::sin(v4_AllN_forM0.Phi() - v4_AllN_forM1.Phi()) * std::sin(phi_MET - v4_AllN_forM1.Phi()) ) ;
    Nu0_pt = pa.Pt();
    Nu1_pt = pb.Pt();
    
    //std::cout<< "----------- " << std::endl;
    //std::cout<< "dphi " << phi_MET - v4_AllN_forM1.Phi() << std::endl;
    //std::cout<< "Nu0_pt " << Nu0_pt << std::endl;
    //std::cout<< "AllN_pt " << v4_AllN_forM0.Pt() << std::endl;

    //Nu1_pt = std::fabs(pt_MET / std::sin(v4_AllN_forM0.Phi() - v4_AllN_forM1.Phi()) * std::sin(phi_MET - v4_AllN_forM0.Phi()) ) ;
    
    dphi0 = (phi_MET - v4_AllN_forM0.Phi())/3.14159265358979323846*180;
    dphi1 = (phi_MET - v4_AllN_forM1.Phi())/3.14159265358979323846*180;
    dphi01 = (v4_AllN_forM0.Phi() - v4_AllN_forM1.Phi())/3.14159265358979323846*180;


    float Nu0Tau_pt = 0;
    float Nu1Tau_pt = 0;

    ROOT::Math::PtEtaPhiMVector paa, pbb;

    DecomposeMomentum(p, v4_tau_forM0, v4_tau_forM1, paa, pbb);

    //Nu0_pt = std::fabs(pt_MET / std::sin(v4_AllN_forM0.Phi() - v4_AllN_forM1.Phi()) * std::sin(phi_MET - v4_AllN_forM1.Phi()) ) ;
    Nu0Tau_pt = paa.Pt();
    Nu1Tau_pt = pbb.Pt();

    mess_dr_01.push_back(Nu0_pt);
    mess_dr_01.push_back(Nu1_pt);
    mess_dr_01.push_back(Nu0Tau_pt);
    mess_dr_01.push_back(Nu1Tau_pt);
    mess_dr_01.push_back(Nu0_pt - v4_AllN_forM0.Pt());
    mess_dr_01.push_back(Nu1_pt - v4_AllN_forM1.Pt());
    mess_dr_01.push_back(Nu0Tau_pt - v4_tau_forM0.Pt());
    mess_dr_01.push_back(Nu1Tau_pt - v4_tau_forM1.Pt());


    //CA
    float _Mass_HTauTau_CA = 0;

    float leading_NuTau_pT_M0 = v4_AllN_forM0.Pt();
    float leading_NuTau_pT_M1 = v4_AllN_forM1.Pt();

    //std::cout<< "matched Nutau0 pt: " << leading_NuTau_pT_M0 << std::endl;
    //std::cout<< "matched Nutau1 pt: " << leading_NuTau_pT_M1 << std::endl;

    float Vis_pT_M0 = GenVisTau_pt[Leading_pT_GenVisTau_Idx[0]];
    float Vis_pT_M1 = GenVisTau_pt[Leading_pT_GenVisTau_Idx[1]];
    
    float x0 = Vis_pT_M0 / (Vis_pT_M0 + leading_NuTau_pT_M0);
    float x1 = Vis_pT_M1 / (Vis_pT_M1 + leading_NuTau_pT_M1);
    
    _Mass_HTauTau_CA = Mass_HTauTau_GenVisTau / std::sqrt(x0 * x1);

    mess_dr_01.push_back(_Mass_HTauTau_CA);


    //MET_to_H
    float x0_METtoNu = Vis_pT_M0 / (Vis_pT_M0 + Nu0_pt);
    float x1_METtoNu = Vis_pT_M1 / (Vis_pT_M1 + Nu1_pt);
    
    float _Mass_HTauTau_CA_METtoNu = Mass_HTauTau_GenVisTau / std::sqrt(x0_METtoNu * x1_METtoNu);

    mess_dr_01.push_back(_Mass_HTauTau_CA_METtoNu);



    float x0_METtoNu_tau = Vis_pT_M0 / (Vis_pT_M0 + Nu0Tau_pt);
    float x1_METtoNu_tau = Vis_pT_M1 / (Vis_pT_M1 + Nu1Tau_pt);
    
    float _Mass_HTauTau_CA_METtoNu_tau = Mass_HTauTau_GenVisTau / std::sqrt(x0_METtoNu_tau * x1_METtoNu_tau);


    //Nu_to_MET
    //float MET_from_Nu01_pt = std::sqrt(v4_AllN_forM0.Pt() * v4_AllN_forM0.Pt() + v4_AllN_forM1.Pt() * v4_AllN_forM1.Pt() 
                                    // + 2 * v4_AllN_forM0.Pt() * v4_AllN_forM1.Pt() * std::cos(v4_AllN_forM0.Phi() - v4_AllN_forM1.Phi()));

    float MET_from_Nu01_pt = (v4_AllN_forM0 + v4_AllN_forM1).Pt();
    float MET_from_Nu01_phi = (v4_AllN_forM0 + v4_AllN_forM1).Phi();


    //float MET_from_Nu01_phi = std::asin(v4_AllN_forM1.Pt() / MET_from_Nu01_pt * std::sin(v4_tau_forM1.Phi() - v4_tau_forM0.Phi())) + v4_tau_forM0.Phi();
    if ( MET_from_Nu01_phi >= M_PI ){
        MET_from_Nu01_phi = MET_from_Nu01_phi - M_PI ;
    }
    else if ( MET_from_Nu01_phi <= -M_PI ){
        MET_from_Nu01_phi = MET_from_Nu01_phi + M_PI ;
    }


    //dMET,need check
    float d_MET_pt = std::sqrt(MET_from_Nu01_pt * MET_from_Nu01_pt + pt_MET * pt_MET 
                            - 2 * pt_MET * MET_from_Nu01_pt * std::cos(MET_from_Nu01_phi - phi_MET));

    float d_MET_phi = std::asin(pt_MET / d_MET_pt * std::sin(phi_MET - MET_from_Nu01_phi)) + MET_from_Nu01_phi;
    if ( d_MET_phi >= M_PI ){
        d_MET_phi = d_MET_phi - M_PI ;
    }
    else if ( d_MET_phi <= -M_PI ){
        d_MET_phi = d_MET_phi + M_PI ;
    }

    //all Nu check
    ROOT::Math::PtEtaPhiMVector v4_AllNu;

    for (int i : GenNutrinoIdx) {

        //int idx = GenNutrinoIdx[i];
        ROOT::Math::PtEtaPhiMVector v4_N(GenPart_pt[i], 
                       GenPart_eta[i], 
                       GenPart_phi[i], 
                       GenPart_mass[i]);
        
        v4_AllNu += v4_N;
    }







    mess_dr_01.push_back(_Mass_HTauTau_CA_METtoNu_tau);

    mess_dr_01.push_back(dphi0);
    mess_dr_01.push_back(dphi1);
    mess_dr_01.push_back(dphi01);

    mess_dr_01.push_back(MET_from_Nu01_pt);
    mess_dr_01.push_back(MET_from_Nu01_phi);

    mess_dr_01.push_back(d_MET_pt);
    mess_dr_01.push_back(d_MET_phi);


    mess_dr_01.push_back(v4_AllNu.Pt());
    mess_dr_01.push_back(v4_AllNu.Phi());

    mess_dr_01.push_back(NuTau_indices_forM0.size());
    mess_dr_01.push_back(NuTau_indices_forM1.size());
    mess_dr_01.push_back(GenNutrinoIdx.size());


    return mess_dr_01;



    
}



float Mass_HTauTau_CA_calc_all_N(float Mass_HTauTau_GenVisTau,
                           cRVecI Leading_pT_GenVisTau_Idx,
                           cRVecI GenVisTau_motherIdx,
                           cRVecI SelGenTauNutrino_motherIdx,
                           cRVecI GenNutrinoIdx,
                           cRVecF GenVisTau_pt,
                           cRVecF GenPart_pt,
                           cRVecF GenPart_eta,
                           cRVecF GenPart_phi,
                           cRVecF GenPart_mass){
    
    float _Mass_HTauTau_CA = 0;
    
    int mother0_idx = GenVisTau_motherIdx[Leading_pT_GenVisTau_Idx[0]];
    int mother1_idx = GenVisTau_motherIdx[Leading_pT_GenVisTau_Idx[1]];
    
    std::vector<int> NuTau_indices_forM0 = findIndicesOfValue(SelGenTauNutrino_motherIdx, mother0_idx);
    std::vector<int> NuTau_indices_forM1 = findIndicesOfValue(SelGenTauNutrino_motherIdx, mother1_idx);

    //std::cout<< "matched Nutau0 number: " << NuTau_indices_forM0.size() << std::endl;
    //std::cout<< "matched Nutau1 number: " << NuTau_indices_forM1.size() << std::endl;
    
    ROOT::Math::PtEtaPhiMVector v4_AllN_forM0;
    ROOT::Math::PtEtaPhiMVector v4_AllN_forM1;

    for (int i : NuTau_indices_forM0) {
        int idx = GenNutrinoIdx[i];
        ROOT::Math::PtEtaPhiMVector v4_N(GenPart_pt[idx], 
                       GenPart_eta[idx], 
                       GenPart_phi[idx], 
                       GenPart_mass[idx]);
        
        v4_AllN_forM0 += v4_N;
    }


    for (int i : NuTau_indices_forM1) {
        int idx = GenNutrinoIdx[i];
        ROOT::Math::PtEtaPhiMVector v4_N(GenPart_pt[idx], 
                       GenPart_eta[idx], 
                       GenPart_phi[idx], 
                       GenPart_mass[idx]);
        
        v4_AllN_forM1 += v4_N;
    }

    
    float leading_NuTau_pT_M0 = v4_AllN_forM0.Pt();
    float leading_NuTau_pT_M1 = v4_AllN_forM1.Pt();

    //std::cout<< "matched Nutau0 pt: " << leading_NuTau_pT_M0 << std::endl;
    //std::cout<< "matched Nutau1 pt: " << leading_NuTau_pT_M1 << std::endl;

    float Vis_pT_M0 = GenVisTau_pt[Leading_pT_GenVisTau_Idx[0]];
    float Vis_pT_M1 = GenVisTau_pt[Leading_pT_GenVisTau_Idx[1]];
    
    float x0 = Vis_pT_M0 / (Vis_pT_M0 + leading_NuTau_pT_M0);
    float x1 = Vis_pT_M1 / (Vis_pT_M1 + leading_NuTau_pT_M1);
    
    _Mass_HTauTau_CA = Mass_HTauTau_GenVisTau / std::sqrt(x0 * x1);
    
    return _Mass_HTauTau_CA;
    
}


float Mass_HTauTau_CA_calc(float Mass_HTauTau_GenVisTau,
                           cRVecI Leading_pT_GenVisTau_Idx,
                           cRVecI GenVisTau_motherIdx,
                           cRVecI SelGenTauNutrino_motherIdx,
                           cRVecI GenNutrinoIdx,
                           cRVecF GenVisTau_pt,
                           cRVecF GenPart_pt){
    
    float _Mass_HTauTau_CA = 0;
    
    int mother0_idx = GenVisTau_motherIdx[Leading_pT_GenVisTau_Idx[0]];
    int mother1_idx = GenVisTau_motherIdx[Leading_pT_GenVisTau_Idx[1]];
    
    std::vector<int> NuTau_indices_forM0 = findIndicesOfValue(SelGenTauNutrino_motherIdx, mother0_idx);
    std::vector<int> NuTau_indices_forM1 = findIndicesOfValue(SelGenTauNutrino_motherIdx, mother1_idx);
    // std::cout<< NuTau_indices_forM0.size() << std::endl;
    RVecF NuTau_forM0_pt;
    for (auto i = 0; i < NuTau_indices_forM0.size(); i++) {

        int idx = GenNutrinoIdx[NuTau_indices_forM0[i]];
        NuTau_forM0_pt.push_back(GenPart_pt[idx]);
    }
    
    RVecF NuTau_forM1_pt;
    for (auto i = 0; i < NuTau_indices_forM1.size(); i++) {
        int idx = GenNutrinoIdx[NuTau_indices_forM1[i]];
        NuTau_forM1_pt.push_back(GenPart_pt[idx]);
    }
    
    auto leadingIter = std::max_element(NuTau_forM0_pt.begin(), NuTau_forM0_pt.end());
    float leading_NuTau_pT_M0 = GenPart_pt[ std::distance(NuTau_forM0_pt.begin(), leadingIter) ];
    // std::cout << leading_NuTau_pT_M0 << std::endl;
    leadingIter = std::max_element(NuTau_forM1_pt.begin(), NuTau_forM1_pt.end());
    float leading_NuTau_pT_M1 = GenPart_pt[ std::distance(NuTau_forM1_pt.begin(), leadingIter) ];
    
    float Vis_pT_M0 = GenVisTau_pt[Leading_pT_GenVisTau_Idx[0]];
    float Vis_pT_M1 = GenVisTau_pt[Leading_pT_GenVisTau_Idx[1]];
    
    float x0 = Vis_pT_M0 / (Vis_pT_M0 + leading_NuTau_pT_M0);
    float x1 = Vis_pT_M1 / (Vis_pT_M1 + leading_NuTau_pT_M1);
    
    _Mass_HTauTau_CA = Mass_HTauTau_GenVisTau / std::sqrt(x0 * x1);
    
    return _Mass_HTauTau_CA;
    // return leading_NuTau_pT_M0;
    
}

//for reco
float min_result(cRVecF vec){
    
    auto leadingIter = *std::min_element(vec.begin(), vec.end());

    return leadingIter;
}

float max_result(cRVecF vec){
    
    auto leadingIter = *std::max_element(vec.begin(), vec.end());

    return leadingIter;
}


int find_tautauJet_from_best_PNetxtt(cRVecF score){
    auto tautauJet_pos = std::distance(score.begin(), std::max_element(score.begin(), score.end()));
    return tautauJet_pos;
}




std::vector<int> find_2tau(cRVecF FatJet_eta, cRVecF FatJet_phi, cRVecF boostedTau_eta, cRVecF boostedTau_phi, int tautauJet_id, cRVecF boostedTau_pt){

    std::vector<int> aaa;

    std::vector<int> dr_taujet_lessthan8;
    std::vector<int> taujet_pt;

    for(auto i = 0; i < boostedTau_eta.size(); i++){
        if (deltaR(FatJet_eta[tautauJet_id],boostedTau_eta[i],FatJet_phi[tautauJet_id],boostedTau_phi[i])<0.8){
            dr_taujet_lessthan8.push_back(i);
            taujet_pt.push_back(boostedTau_pt[i]);
        }
    }
    std::vector<int> indices(dr_taujet_lessthan8.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        indices[i] = i;
    }
    std::sort(indices.begin(), indices.end(), [&taujet_pt](int a, int b) {
        return taujet_pt[a] > taujet_pt[b];
    });

    std::vector<int> sorteddr_taujet_lessthan8(dr_taujet_lessthan8.size());
    std::vector<int> sortedtaujet_pt(taujet_pt.size());

    for (size_t i = 0; i < indices.size(); ++i) {
        sorteddr_taujet_lessthan8[i] = dr_taujet_lessthan8[indices[i]];
        sortedtaujet_pt[i] = taujet_pt[indices[i]];
    }

    dr_taujet_lessthan8 = sorteddr_taujet_lessthan8;
    taujet_pt = sortedtaujet_pt;


    aaa.push_back(dr_taujet_lessthan8.size());
    aaa.insert(aaa.end(), dr_taujet_lessthan8.begin(), dr_taujet_lessthan8.end());
    return aaa;
}

float M_tautau_CA(cRVecF boostedTau_eta, cRVecF boostedTau_phi, 
                    cRVecF boostedTau_mass, cRVecF boostedTau_pt, 
                    int tau0_id, int tau1_id,
                    float GenMET_pt,float GenMET_phi,
                    float mass_vistautau){
    float pt_MET = GenMET_pt;
    float phi_MET = GenMET_phi;


    ROOT::Math::PtEtaPhiMVector v4_AllN_forM0(boostedTau_pt[tau0_id], 
                       boostedTau_eta[tau0_id], 
                       boostedTau_phi[tau0_id], 
                       boostedTau_mass[tau0_id]);
    
    ROOT::Math::PtEtaPhiMVector v4_AllN_forM1(boostedTau_pt[tau1_id], 
                       boostedTau_eta[tau1_id], 
                       boostedTau_phi[tau1_id], 
                       boostedTau_mass[tau1_id]);

    ROOT::Math::PtEtaPhiMVector p_MET(pt_MET, 0, phi_MET, pt_MET);
    ROOT::Math::PtEtaPhiMVector p_Nu0, p_Nu1;

    DecomposeMomentum(p_MET, v4_AllN_forM0, v4_AllN_forM1, p_Nu0, p_Nu1);

    float Nu0Tau_pt = p_Nu0.Pt();
    float Nu1Tau_pt = p_Nu1.Pt();

    float x0_METtoNu_tau = boostedTau_pt[tau0_id] / (boostedTau_pt[tau0_id] + Nu0Tau_pt);
    float x1_METtoNu_tau = boostedTau_pt[tau1_id] / (boostedTau_pt[tau1_id] + Nu1Tau_pt);
    
    float M_tautau_CA = mass_vistautau / std::sqrt(x0_METtoNu_tau * x1_METtoNu_tau);

    return M_tautau_CA;

    }

int veto_e_muon(cRVecF Muon_pt, cRVecF Muon_eta, cRVecF Muon_dxy, cRVecF Muon_dz, float pt_cut, float eta_cut, float dxy_cut, float dz_cut){
    
    int e_muon_number = 0;

    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if((Muon_pt[i] > pt_cut) && (std::abs(Muon_eta[i]) < eta_cut) && (std::abs(Muon_dxy[i]) < dxy_cut) && (std::abs(Muon_dz[i]) < dz_cut)){
            e_muon_number += 1;
        }
    }

    return e_muon_number;
}

