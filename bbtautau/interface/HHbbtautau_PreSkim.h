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


RVecI GenTauNutrino_Selector(cRVecI GenPart_pdgId){

    RVecI GenNutrinoIdx;
    
    int idx = -99;
    for(auto i = 0; i < GenPart_pdgId.size(); i++) {
        if ( (abs(GenPart_pdgId[i]) == 16) ) {
            idx = i;
            GenNutrinoIdx.push_back(idx);
        }
    }

    return GenNutrinoIdx; // GenTauIdx(idx_nvtau1, idx_nvtau2, ...)
}


RVecI SelGenTauNutrino_Mother(cRVecI GenTauNutrino_idx , 
                              cRVecI GenPart_pdgId,
                              cRVecF GenPart_genPartIdxMother, 
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
                       cRVecF GenPart_genPartIdxMother, 
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

    int _nMatch = 0;

    for (auto i = 0; i < SelGenTauNutrino_motherIdx.size(); i++) {
        if (isNumberInList(SelGenTau_Mother_idx, SelGenTauNutrino_motherIdx[i])) {
            _nMatch++;
        }
    }
    
    return _nMatch;
        
}

int GenTauNutrino_GenVisTau_Match(cRVecI SelGenTauNutrino_motherIdx, 
                                  cRVecI GenVisTau_motherIdx){

    int _nMatch = 0;

    for (auto i = 0; i < SelGenTauNutrino_motherIdx.size(); i++) {
        if (isNumberInList(GenVisTau_motherIdx, SelGenTauNutrino_motherIdx[i])) {
            _nMatch++;
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



float Mass_HTauTau_CA_calc(float Mass_HTauTau_GenVisTau,
                           cRVecI Leading_pT_GenVisTau_Idx,
                           cRVecI GenVisTau_motherIdx,
                           cRVecI SelGenTauNutrino_motherIdx,
                           cRVecF GenVisTau_pt,
                           cRVecF GenPart_pt){
    
    float _Mass_HTauTau_CA = 0;
    
    int mother0_idx = GenVisTau_motherIdx[Leading_pT_GenVisTau_Idx[0]];
    int mother1_idx = GenVisTau_motherIdx[Leading_pT_GenVisTau_Idx[1]];
    
    std::vector<int> NuTau_indices_forM0 = findIndicesOfValue(SelGenTauNutrino_motherIdx, mother0_idx);
    std::vector<int> NuTau_indices_forM1 = findIndicesOfValue(SelGenTauNutrino_motherIdx, mother1_idx);
    
    RVecF NuTau_forM0_pt;
    for (auto i = 0; i < NuTau_indices_forM0.size(); i++) {
        NuTau_forM0_pt.push_back(GenPart_pt[NuTau_indices_forM0[i]]);
    }
    
    RVecF NuTau_forM1_pt;
    for (auto i = 0; i < NuTau_indices_forM1.size(); i++) {
        NuTau_forM1_pt.push_back(GenPart_pt[NuTau_indices_forM1[i]]);
    }
    
    auto leadingIter = std::max_element(NuTau_forM0_pt.begin(), NuTau_forM0_pt.end());
    float leading_NuTau_pT_M0 = GenPart_pt[ std::distance(NuTau_forM0_pt.begin(), leadingIter) ];
    
    leadingIter = std::max_element(NuTau_forM1_pt.begin(), NuTau_forM1_pt.end());
    float leading_NuTau_pT_M1 = GenPart_pt[ std::distance(NuTau_forM1_pt.begin(), leadingIter) ];
    
    float Vis_pT_M0 = GenVisTau_pt[Leading_pT_GenVisTau_Idx[0]];
    float Vis_pT_M1 = GenVisTau_pt[Leading_pT_GenVisTau_Idx[1]];
    
    float x0 = Vis_pT_M0 / (Vis_pT_M0 + leading_NuTau_pT_M0);
    float x1 = Vis_pT_M1 / (Vis_pT_M1 + leading_NuTau_pT_M1);
    
    _Mass_HTauTau_CA = Mass_HTauTau_GenVisTau / std::sqrt(x0 * x1);
    
    return _Mass_HTauTau_CA;
    
}
