import sys
import os
import math
import argparse
import ROOT as R
import numpy as np
from array import *
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--inputDirectory",          dest="ind",   default="test_0908_catfirst_bdt_30bins/",   help="input directory" )
parser.add_argument("-o", "--output",          dest="out",   default="",   help="output name" )
parser.add_argument("-y", "--year",          dest="year",   default="2018",   help="output name" )
args = parser.parse_args()

# Manipulate histograms
# If the first and last bin is too large, make it short 
# If the first x bins and last y bins are empty, merge it to adjacent bins

if args.ind != 0:
    path = str(args.ind)+'/'

import os
data_outDir = args.out
if data_outDir is not '' :
    if not os.path.isdir( data_outDir ) :
        os.makedirs( data_outDir )

syst = {}
addtional = ''


data_obs = {}#"data_obs"}

signal_procs = ["VHH_CV_0p5_C2V_1_kl_1_hbbhbb", "VHH_CV_1_C2V_0_kl_1_hbbhbb", "VHH_CV_1_C2V_1_kl_1_hbbhbb", "VHH_CV_1_C2V_1_kl_2_hbbhbb", "VHH_CV_1_C2V_2_kl_1_hbbhbb", "VHH_CV_1p5_C2V_1_kl_1_hbbhbb","VHH_CV_1_C2V_1_kl_0_hbbhbb",'VHH_CV_1_C2V_1_kl_20_hbbhbb',
"WHH_CV_0p5_C2W_1_kl_1_hbbhbb", "WHH_CV_1_C2W_0_kl_1_hbbhbb", "WHH_CV_1_C2W_1_kl_1_hbbhbb", "WHH_CV_1_C2W_1_kl_2_hbbhbb", "WHH_CV_1_C2W_2_kl_1_hbbhbb", "WHH_CV_1p5_C2W_1_kl_1_hbbhbb","WHH_CV_1_C2W_1_kl_0_hbbhbb",'WHH_CV_1_C2W_1_kl_20_hbbhbb',
"ZHH_CV_0p5_C2Z_1_kl_1_hbbhbb", "ZHH_CV_1_C2Z_0_kl_1_hbbhbb", "ZHH_CV_1_C2Z_1_kl_1_hbbhbb", "ZHH_CV_1_C2Z_1_kl_2_hbbhbb", "ZHH_CV_1_C2Z_2_kl_1_hbbhbb", "ZHH_CV_1p5_C2Z_1_kl_1_hbbhbb","ZHH_CV_1_C2Z_1_kl_0_hbbhbb",'ZHH_CV_1_C2Z_1_kl_20_hbbhbb',
]


plots_proce = ["TT_FailReweight","TTB_FailReweight","VHH_CV_1_C2V_1_kl_1_hbbhbb","VHH_CV_1_C2V_1_kl_2_hbbhbb","VHH_CV_1_C2V_1_kl_0_hbbhbb"]
plots_proce = ["VHH_CV_1_C2V_1_kl_1_hbbhbb","VHH_CV_1_C2V_1_kl_2_hbbhbb","VHH_CV_1_C2V_1_kl_0_hbbhbb"]

VHHfiles = [i for i in os.listdir(path) if i.endswith(".root")]

procs = ["VHH_CV_0p5_C2V_1_kl_1_hbbhbb", "VHH_CV_1_C2V_0_kl_1_hbbhbb", "VHH_CV_1_C2V_1_kl_1_hbbhbb", "VHH_CV_1_C2V_1_kl_2_hbbhbb", "VHH_CV_1_C2V_2_kl_1_hbbhbb", "VHH_CV_1p5_C2V_1_kl_1_hbbhbb","TT","TTB","data_obs", 'VHH_CV_1_C2V_1_kl_0_hbbhbb','VHH_CV_1_C2V_1_kl_20_hbbhbb','ST','ttH_hbb','TTV','Znunu','DY', "WHH_CV_0p5_C2W_1_kl_1_hbbhbb", "WHH_CV_1_C2W_0_kl_1_hbbhbb", "WHH_CV_1_C2W_1_kl_1_hbbhbb", "WHH_CV_1_C2W_1_kl_2_hbbhbb", "WHH_CV_1_C2W_2_kl_1_hbbhbb", "WHH_CV_1p5_C2W_1_kl_1_hbbhbb","WHH_CV_1_C2W_1_kl_0_hbbhbb",'WHH_CV_1_C2W_1_kl_20_hbbhbb',
"ZHH_CV_0p5_C2Z_1_kl_1_hbbhbb", "ZHH_CV_1_C2Z_0_kl_1_hbbhbb", "ZHH_CV_1_C2Z_1_kl_1_hbbhbb", "ZHH_CV_1_C2Z_1_kl_2_hbbhbb", "ZHH_CV_1_C2Z_2_kl_1_hbbhbb", "ZHH_CV_1p5_C2Z_1_kl_1_hbbhbb","ZHH_CV_1_C2Z_1_kl_0_hbbhbb",'ZHH_CV_1_C2Z_1_kl_20_hbbhbb',
]

plots_proce = ["TT","TTB","VHH_CV_1_C2V_1_kl_1_hbbhbb","VHH_CV_1_C2V_1_kl_2_hbbhbb","VHH_CV_1_C2V_1_kl_0_hbbhbb",'ST']
bkg_procs = {"TT","TTB"}

colorb = [R.kRed,R.kBlue,R.kGreen,R.kBlack,R.kCyan,R.kOrange,R.kViolet,R.kGray,R.kPink]


def replace_hists():
    friendpath='test_221005_collectpileup_btag_LHE/'
    for ifile in VHHfiles:
        print(ifile)
        hsample = {}
        hout={}
        fin = R.TFile.Open('{0}'.format(path+ifile),'READ')
        alist=fin.GetListOfKeys()
        alistname=[]
        for i in range(alist.GetEntries()):
            alistname.append(alist.At(i).GetName())

        friend = R.TFile.Open('{0}'.format(friendpath+ifile),'READ')
        alistfriend=friend.GetListOfKeys()
        alistnamefirend=[]
        for i in range(alistfriend.GetEntries()):
            alistnamefirend.append(alistfriend.At(i).GetName())

        for iprocs in alistname:
            if ( "CMS_pileup" in iprocs  or 'CMS_btag_LF_2016_2017_2018' in iprocs or 'CMS_btag_cferr1_2016_2017_2018' in iprocs or 'CMS_btag_cferr2_2016_2017_2018' in iprocs) and 'true' not in iprocs:
            #if False:
                print('gonna replace ',iprocs,' in ', ifile)
                friend.cd()
                hout[iprocs] = friend.Get(iprocs).Clone()
            else:
                fin.cd()
                hout[iprocs] = fin.Get(iprocs).Clone()

        for iprocs in alistnamefirend:
            if iprocs not in alistname:
                print("I don't see ", iprocs, "add it from friend?")
                hout[iprocs] = friend.Get(iprocs).Clone()

        fout = R.TFile.Open('{0}/{1}{2}'.format(args.out,addtional,ifile),'RECREATE')
        for iprocs in hout.keys():
            hout[iprocs].Write()
        fout.Close()

#replace_hists()
#exit()

def merge_signals():
    WHH=["WHH_CV_0p5_C2W_1_kl_1_hbbhbb", "WHH_CV_1_C2W_0_kl_1_hbbhbb", "WHH_CV_1_C2W_1_kl_1_hbbhbb", "WHH_CV_1_C2W_1_kl_2_hbbhbb", "WHH_CV_1_C2W_2_kl_1_hbbhbb", "WHH_CV_1p5_C2W_1_kl_1_hbbhbb","WHH_CV_1_C2W_1_kl_0_hbbhbb",'WHH_CV_1_C2W_1_kl_20_hbbhbb']
    ZHH=["ZHH_CV_0p5_C2Z_1_kl_1_hbbhbb", "ZHH_CV_1_C2Z_0_kl_1_hbbhbb", "ZHH_CV_1_C2Z_1_kl_1_hbbhbb", "ZHH_CV_1_C2Z_1_kl_2_hbbhbb", "ZHH_CV_1_C2Z_2_kl_1_hbbhbb", "ZHH_CV_1p5_C2Z_1_kl_1_hbbhbb","ZHH_CV_1_C2Z_1_kl_0_hbbhbb",'ZHH_CV_1_C2Z_1_kl_20_hbbhbb',]
    for ifile in VHHfiles:
        print(ifile)
        hsample = {}
        hout={}
        fin = R.TFile.Open('{0}'.format(path+ifile),'READ')
        hist_template = fin.Get('ST').Clone()
        hist_template.Reset()
        alist=fin.GetListOfKeys()
        alistname=[]
        for i in range(alist.GetEntries()):
            alistname.append(alist.At(i).GetName())        
        fin.cd()
        for iprocs in alistname:
            if 'HH' not in iprocs: hout[iprocs] = fin.Get(iprocs).Clone()
            else:
                if iprocs.replace("WHH_CV","VHH_CV").replace("C2W","C2V").replace("ZHH_CV","VHH_CV").replace("C2Z","C2V") in hout.keys(): continue
                for isig in WHH:
                    if isig in iprocs: 
                        hout[iprocs.replace("WHH_CV","VHH_CV").replace("C2W","C2V")] = fin.Get(iprocs).Clone()
                        if iprocs.replace("WHH_CV","ZHH_CV").replace("C2W","C2Z") in alistname: hout[iprocs.replace("WHH_CV","VHH_CV").replace("C2W","C2V")].Add(fin.Get(iprocs.replace("WHH_CV","ZHH_CV").replace("C2W","C2Z")))
                        elif isig.replace("WHH_CV","ZHH_CV").replace("C2W","C2Z") in alistname: hout[iprocs.replace("WHH_CV","VHH_CV").replace("C2W","C2V")].Add(fin.Get(isig.replace("WHH_CV","ZHH_CV").replace("C2W","C2Z")))
                        hout[iprocs.replace("WHH_CV","VHH_CV").replace("C2W","C2V")].SetName(iprocs.replace("WHH_CV","VHH_CV").replace("C2W","C2V"))
                for isig in ZHH:
                    if isig in iprocs:
                        hout[iprocs.replace("ZHH_CV","VHH_CV").replace("C2Z","C2V")] = fin.Get(iprocs).Clone()
                        if iprocs.replace("ZHH_CV","WHH_CV").replace("C2Z","C2Z") in alistname: hout[iprocs.replace("ZHH_CV","VHH_CV").replace("C2W","C2V")].Add(fin.Get(iprocs.replace("ZHH_CV","WHH_CV").replace("C2Z","C2W")))
                        elif isig.replace("ZHH_CV","WHH_CV").replace("C2Z","C2W") in alistname: hout[iprocs.replace("ZHH_CV","VHH_CV").replace("C2Z","C2V")].Add(fin.Get(isig.replace("ZHH_CV","WHH_CV").replace("C2Z","C2W")))
                        hout[iprocs.replace("ZHH_CV","VHH_CV").replace("C2Z","C2V")].SetName(iprocs.replace("ZHH_CV","VHH_CV").replace("C2Z","C2V"))

        
        # merge PDF uncertainties
        # TT quadratic sum of all residuals 
        # Others RMS
        mergePDF=False
        if mergePDF:
            strategy=2
            binning=[]
            for i in range(1,hist_template.GetNbinsX()+2):
                binning.append(hist_template.GetBinLowEdge(i))

            processes_with_PDF=[]
            for iprocs in hout.keys():
                if "CMS_vhh_PDF" in iprocs and iprocs.split("_CMS_vhh_PDF")[0] not in processes_with_PDF: processes_with_PDF.append(iprocs.split("_CMS_vhh_PDF")[0])
                else:    continue
            print('need to worry about the PDF of ', processes_with_PDF)
            PDFname={
                'TT':'_CMS_vhh_PDF_tt_',
                'ST':'_CMS_vhh_PDF_st',
                'TTB':'_CMS_vhh_PDF_ttbb_',
                'TTV':'_CMS_vhh_PDF_ttv_',
                'ttH':'_CMS_vhh_PDF_ttH_',
            }        
            for iprocess in signal_procs:
                PDFname[iprocess]='_CMS_vhh_PDF_'
            for iprocess in processes_with_PDF:
                hout[iprocess+'_CMS_vhh_PDFUp'] = R.TH1F(iprocess+'_CMS_vhh_PDFUp',iprocess+'_CMS_vhh_PDFUp',len(binning)-1,array('d',binning))
                hout[iprocess+'_CMS_vhh_PDFDown'] = R.TH1F(iprocess+'_CMS_vhh_PDFDown',iprocess+'_CMS_vhh_PDFDown',len(binning)-1,array('d',binning))
                for i in range(len(binning)-1):
                    pdfvaries=[]
                    for j in range(50):
                        pdfvaries.append(hout[iprocess+PDFname[iprocess]+str(j)+'Up'].GetBinContent(i+1))
                        pdfvaries.append(hout[iprocess+PDFname[iprocess]+str(j)+'Down'].GetBinContent(i+1))
                    if strategy==1:
                        pdfvariesmax=0
                        for j in range(len(pdfvaries)):
                            shift = abs(pdfvaries[j] - hout[iprocess].GetBinContent(i+1))
                            if shift>pdfvariesmax: pdfvariesmax=shift
                        hout[iprocess+'_CMS_vhh_PDFUp'].SetBinContent(i+1, hout[iprocess].GetBinContent(i+1) +pdfvariesmax )
                        hout[iprocess+'_CMS_vhh_PDFDown'].SetBinContent(i+1, hout[iprocess].GetBinContent(i+1)-pdfvariesmax)
                        if hout[iprocess].GetBinContent(i+1)!=0:
                            hout[iprocess+'_CMS_vhh_PDFUp'].SetBinError(i+1, hout[iprocess].GetBinError(i+1)*(pdfvariesmax+hout[iprocess].GetBinContent(i+1))/hout[iprocess].GetBinContent(i+1) )
                            hout[iprocess+'_CMS_vhh_PDFDown'].SetBinError(i+1, hout[iprocess].GetBinError(i+1)*(-pdfvariesmax+hout[iprocess].GetBinContent(i+1))/hout[iprocess].GetBinContent(i+1) )
                        else:
                            hout[iprocess+'_CMS_vhh_PDFUp'].SetBinError(i+1, hout[iprocess].GetBinError(i+1) )
                            hout[iprocess+'_CMS_vhh_PDFDown'].SetBinError(i+1, hout[iprocess].GetBinError(i+1) )
                    else:
                        if iprocess!="TTB": 
                            pdfvariestot=0
                            for j in range(len(pdfvaries)):
                                shift = abs(pdfvaries[j] - hout[iprocess].GetBinContent(i+1))
                                pdfvariestot += shift**2
                            hout[iprocess+'_CMS_vhh_PDFUp'].SetBinContent(i+1, np.sqrt(pdfvariestot)+hout[iprocess].GetBinContent(i+1) )
                            hout[iprocess+'_CMS_vhh_PDFDown'].SetBinContent(i+1, hout[iprocess].GetBinContent(i+1)-np.sqrt(pdfvariestot) )
                            if hout[iprocess].GetBinContent(i+1)!=0:
                                hout[iprocess+'_CMS_vhh_PDFUp'].SetBinError(i+1, hout[iprocess].GetBinError(i+1)*(np.sqrt(pdfvariestot)+hout[iprocess].GetBinContent(i+1))/hout[iprocess].GetBinContent(i+1) )
                                hout[iprocess+'_CMS_vhh_PDFDown'].SetBinError(i+1, hout[iprocess].GetBinError(i+1)*(-np.sqrt(pdfvariestot)+hout[iprocess].GetBinContent(i+1))/hout[iprocess].GetBinContent(i+1) )
                            else:
                                hout[iprocess+'_CMS_vhh_PDFUp'].SetBinError(i+1, hout[iprocess].GetBinError(i+1) )
                                hout[iprocess+'_CMS_vhh_PDFDown'].SetBinError(i+1, hout[iprocess].GetBinError(i+1) )
                        else:
                            hout[iprocess+'_CMS_vhh_PDFUp'].SetBinContent(i+1, np.std(pdfvaries)+hout[iprocess].GetBinContent(i+1) )
                            hout[iprocess+'_CMS_vhh_PDFDown'].SetBinContent(i+1, hout[iprocess].GetBinContent(i+1)-np.std(pdfvaries) )
                            if hout[iprocess].GetBinContent(i+1):
                                hout[iprocess+'_CMS_vhh_PDFUp'].SetBinError(i+1, hout[iprocess].GetBinError(i+1)*(np.std(pdfvaries)+hout[iprocess].GetBinContent(i+1))/hout[iprocess].GetBinContent(i+1) )
                                hout[iprocess+'_CMS_vhh_PDFDown'].SetBinError(i+1, hout[iprocess].GetBinError(i+1)*(-np.std(pdfvaries)+hout[iprocess].GetBinContent(i+1))/hout[iprocess].GetBinContent(i+1) )
                            else:
                                hout[iprocess+'_CMS_vhh_PDFUp'].SetBinError(i+1, hout[iprocess].GetBinError(i+1) )
                                hout[iprocess+'_CMS_vhh_PDFDown'].SetBinError(i+1, hout[iprocess].GetBinError(i+1) )

            processes_with_scaleMuRMuF=[]
            for iprocs in hout.keys():
                if "CMS_vhh_scaleMuRMuF" in iprocs and iprocs.split("_CMS_vhh_scaleMuRMuF")[0] not in processes_with_scaleMuRMuF: processes_with_scaleMuRMuF.append(iprocs.split("_CMS_vhh_scaleMuRMuF")[0])
                else:    continue
            print('need to worry about the scaleMuRMuF of ', processes_with_scaleMuRMuF)
            scaleMuRMuFname={
                'TT':'_CMS_vhh_scaleMuRMuF_tt_',
                'ST':'_CMS_vhh_scaleMuRMuF_st_',
                'TTB':'_CMS_vhh_scaleMuRMuF_ttbb_',
                'TTV':'_CMS_vhh_scaleMuRMuF_ttv_',
                'ttH':'_CMS_vhh_scaleMuRMuF_ttH_',
            }
            strategy=1
            for iprocess in signal_procs:
                scaleMuRMuFname[iprocess]='_CMS_vhh_scaleMuRMuF_'
            for iprocess in processes_with_scaleMuRMuF:
                hout[iprocess+'_CMS_vhh_scaleMuRMuFUp'] = R.TH1F(iprocess+'_CMS_vhh_scaleMuRMuFUp',iprocess+'_CMS_vhh_scaleMuRMuFUp',len(binning)-1,array('d',binning))
                hout[iprocess+'_CMS_vhh_scaleMuRMuFDown'] = R.TH1F(iprocess+'_CMS_vhh_scaleMuRMuFDown',iprocess+'_CMS_vhh_scaleMuRMuFDown',len(binning)-1,array('d',binning))
                for i in range(len(binning)-1):
                    pdfvaries=[]
                    for j in range(4):
                        pdfvaries.append(hout[iprocess+scaleMuRMuFname[iprocess]+str(j)+'Up'].GetBinContent(i+1))
                        pdfvaries.append(hout[iprocess+scaleMuRMuFname[iprocess]+str(j)+'Down'].GetBinContent(i+1))
                    if strategy==1:
                        hout[iprocess+'_CMS_vhh_scaleMuRMuFUp'].SetBinContent(i+1, max(pdfvaries) )
                        hout[iprocess+'_CMS_vhh_scaleMuRMuFDown'].SetBinContent(i+1,min(pdfvaries))
                        if hout[iprocess].GetBinContent(i+1)!=0:
                            hout[iprocess+'_CMS_vhh_scaleMuRMuFUp'].SetBinError(i+1, hout[iprocess].GetBinError(i+1)*max(pdfvaries) /hout[iprocess].GetBinContent(i+1) )
                            hout[iprocess+'_CMS_vhh_scaleMuRMuFDown'].SetBinError(i+1, hout[iprocess].GetBinError(i+1)*min(pdfvaries) /hout[iprocess].GetBinContent(i+1) )

        fout = R.TFile.Open('{0}/{1}{2}'.format(args.out,addtional,ifile),'RECREATE')
        for iprocs in hout.keys():
            if 'CMS_vhh4b_scaleMuRMuF_ttH' in iprocs: hout[iprocs].SetName(iprocs.replace('ttH_CMS',"ttH_hbb_CMS"))
            elif 'ttH' in iprocs and 'ttH_hbb' not in iprocs:  hout[iprocs].SetName(iprocs.replace("ttH",'ttH_hbb'))
            hout[iprocs].Write()
        fout.Close()

def merge_years():
    cards=[]
    for ifile in VHHfiles:
        newfile= ifile.replace("2016","run2").replace("2017","run2").replace("2018","run2")
        if newfile in cards or 'Zll' not in ifile: 
            # already processed
            continue
        else: cards.append(newfile)
        print(ifile)
        # start process newfile
        hsample = {}
        hout={}
        filelist={}
        finlist={}
        alistname={}
        for iyear in ['2016','2017','2018']:
            filelist[iyear]=newfile.replace("run2",iyear)
            finlist[iyear] = R.TFile.Open('{0}'.format(path+filelist[iyear]),'READ')
            alist=finlist[iyear].GetListOfKeys()
            alistname[iyear]=[]
            for i in range(alist.GetEntries()):
                alistname[iyear].append(alist.At(i).GetName())

        all_process=[]
        for iprocs in alistname['2016']+alistname['2017']+alistname['2018']:
            if 'true' in iprocs: continue
            if iprocs in all_process: continue
            #if 'TT' not in iprocs: continue
            #if 'TTB' in iprocs or 'TTV' in iprocs: continue
            print(" ")
            if iprocs in procs:
                print('this is norminal', iprocs)
                print(finlist['2016'].Get(iprocs).Integral(), finlist['2017'].Get(iprocs).Integral(),finlist['2018'].Get(iprocs).Integral(),)
                finlist['2016'].cd()
                hout[iprocs] = finlist['2016'].Get(iprocs).Clone()
                finlist['2017'].cd()
                hout[iprocs].Add(finlist['2017'].Get(iprocs))
                finlist['2018'].cd()
                hout[iprocs].Add(finlist['2018'].Get(iprocs))
                print(hout[iprocs].Integral())
                all_process.append(iprocs)
            elif '2016_2017_2018' in iprocs: 
                print('this should be correlated', iprocs)
                print(finlist['2016'].Get(iprocs).Integral(), finlist['2017'].Get(iprocs).Integral(),finlist['2018'].Get(iprocs).Integral(),)
                finlist['2016'].cd()
                hout[iprocs] = finlist['2016'].Get(iprocs).Clone()
                finlist['2017'].cd()
                hout[iprocs].Add(finlist['2017'].Get(iprocs))
                finlist['2018'].cd()
                hout[iprocs].Add(finlist['2018'].Get(iprocs))
                print(hout[iprocs].Integral())
                all_process.append(iprocs)
            elif '2016' in iprocs:
                print(iprocs, 'is 2016 specific, will use ', iprocs.split('_CMS')[0], ' for 17/18')
                print(finlist['2016'].Get(iprocs).Integral(), finlist['2017'].Get(iprocs.split('_CMS')[0]).Integral(),finlist['2018'].Get(iprocs.split('_CMS')[0]).Integral(),)
                finlist['2016'].cd()
                hout[iprocs] = finlist['2016'].Get(iprocs).Clone()
                finlist['2017'].cd()
                hout[iprocs].Add(finlist['2017'].Get(iprocs.split('_CMS')[0]))
                finlist['2018'].cd()
                hout[iprocs].Add(finlist['2018'].Get(iprocs.split('_CMS')[0]))
                all_process.append(iprocs)
            elif '2017' in iprocs:
                print(iprocs, 'is 2017 specific, will use ', iprocs.split('_CMS')[0], ' for 16/18')
                finlist['2017'].cd()
                hout[iprocs] = finlist['2017'].Get(iprocs).Clone()
                finlist['2016'].cd()
                hout[iprocs].Add(finlist['2016'].Get(iprocs.split('_CMS')[0]))
                finlist['2018'].cd()
                hout[iprocs].Add(finlist['2018'].Get(iprocs.split('_CMS')[0]))
                all_process.append(iprocs)
            elif '2018' in iprocs:
                print(iprocs, 'is 2018 specific, will use ', iprocs.split('_CMS')[0], ' for 16/17')
                finlist['2018'].cd()
                hout[iprocs] = finlist['2018'].Get(iprocs).Clone()
                finlist['2016'].cd()
                hout[iprocs].Add(finlist['2016'].Get(iprocs.split('_CMS')[0]))
                finlist['2017'].cd()
                hout[iprocs].Add(finlist['2017'].Get(iprocs.split('_CMS')[0]))
                all_process.append(iprocs)
            elif 'CMS' in iprocs:
                print('should also correlate ', iprocs)
                print(finlist['2016'].Get(iprocs).Integral(), finlist['2017'].Get(iprocs).Integral(),finlist['2018'].Get(iprocs).Integral(),)
                finlist['2016'].cd()
                hout[iprocs] = finlist['2016'].Get(iprocs).Clone()
                finlist['2017'].cd()
                hout[iprocs].Add(finlist['2017'].Get(iprocs))
                finlist['2018'].cd()
                hout[iprocs].Add(finlist['2018'].Get(iprocs))
                all_process.append(iprocs)
            else:
                print('!!', iprocs)
        fout = R.TFile.Open('{0}/{1}{2}'.format(args.out,addtional,newfile),'RECREATE')
        #if '3b' in newfile:
        #    for isys in ['Up','Down']:
        #        hout['DY_CMS_vhh4b_DYZptCor'+isys] = hout['DY_CMS_vhh4b_DYZptCor_3b'+isys].Clone()
        #        hout['DY_CMS_vhh4b_DYZptCor'+isys].SetName('DY_CMS_vhh4b_DYZptCor'+isys)
        #        print(hout['DY_CMS_vhh4b_DYZptCor_3b'+isys].Integral(), hout['DY_CMS_vhh4b_DYZptCor'+isys].Integral())
        #elif '4b' in newfile:
        #    for isys in ['Up','Down']:
        #        hout['DY_CMS_vhh4b_DYZptCor'+isys] = hout['DY_CMS_vhh4b_DYZptCor_4b'+isys].Clone()
        #        hout['DY_CMS_vhh4b_DYZptCor'+isys].SetName('DY_CMS_vhh4b_DYZptCor'+isys)
        #        print(hout['DY_CMS_vhh4b_DYZptCor_4b'+isys].Integral(), hout['DY_CMS_vhh4b_DYZptCor'+isys].Integral())
        #else:
        #    print('no DYZptCor?')
        #    print(hout.keys())
        #    exit()
        for iprocs in hout.keys():
            hout[iprocs].Write()
        fout.Close()


def merge_ZllTTCR_3b4b():
    cards=[]
    VHHfiles=['sum_hists_run2_R_Zll_3b_TTCR_.root', 'sum_hists_run2_R_Zll_4b_TTCR_.root']
    for ifile in VHHfiles:
        newfile= ifile.replace("3b","mergeb").replace("4b","mergeb")
        if newfile in cards or 'Zll' not in ifile:
            continue
        else: cards.append(newfile)
        print(ifile)
        # start process newfile
        hsample = {}
        hout={}
        filelist={}
        finlist={}
        alistname={}
        for iyear in ['3b','4b']:
            filelist[iyear]=newfile.replace("mergeb",iyear)
            finlist[iyear] = R.TFile.Open('{0}'.format(path+filelist[iyear]),'READ')
            alist=finlist[iyear].GetListOfKeys()
            alistname[iyear]=[]
            for i in range(alist.GetEntries()):
                alistname[iyear].append(alist.At(i).GetName())

        all_process=[]
        for iprocs in alistname['3b']+alistname['4b']:
            if 'true' in iprocs: continue
            if iprocs in all_process: continue
            #if "TT_CMS_eff_ee_2018Up" == iprocs: print(iprocs)
            #if "TT" not in iprocs: continue
            if 'bkgrwt' in iprocs and 'R_3b' in iprocs:
                print('this is rwt', iprocs)
                finlist['3b'].cd()
                hout[iprocs] = finlist['3b'].Get(iprocs).Clone()
                finlist['4b'].cd()
                hout[iprocs].Add(finlist['4b'].Get(iprocs.split('_CMS')[0]))
                all_process.append(iprocs)
            elif 'bkgrwt' in iprocs and 'R_4b' in iprocs:
                print('this is rwt', iprocs)
                finlist['3b'].cd()
                hout[iprocs] = finlist['4b'].Get(iprocs).Clone()
                finlist['4b'].cd()
                hout[iprocs].Add(finlist['3b'].Get(iprocs.split('_CMS')[0]))
                all_process.append(iprocs)
            else:
                #print(iprocs, finlist['3b'].Get(iprocs).Integral(), finlist['4b'].Get(iprocs).Integral())
                finlist['3b'].cd()
                hout[iprocs] = finlist['3b'].Get(iprocs).Clone()
                finlist['4b'].cd()
                hout[iprocs].Add(finlist['4b'].Get(iprocs))
                #print(hout[iprocs].Integral())
                all_process.append(iprocs)
        fout = R.TFile.Open('{0}/{1}{2}'.format(args.out,addtional,newfile),'RECREATE')
        for iprocs in hout.keys():
            if 'TT_' in iprocs: print(iprocs, hout[iprocs].Integral())
            hout[iprocs].Write()
        fout.Close()


def merge_emptyBins():
    for ifile in VHHfiles:
        print(ifile)
        hsample = {}
        hout={}
        fin = R.TFile.Open('{0}'.format(path+ifile),'READ')
        hist_template = fin.Get('ST').Clone()
        hist_template.Reset()
        binning=[]
        for i in range(1,hist_template.GetNbinsX()+2):
            binning.append(hist_template.GetBinLowEdge(i))
        EMPTYN=0
        EMPTYX=[]
        bkg_procs = {"TTB","TT"}
        alist=fin.GetListOfKeys()
        alistname=[]
        for i in range(alist.GetEntries()):
            alistname.append(alist.At(i).GetName())
        fin.cd()
        hist_template = fin.Get("TT").Clone()
        hist_template.Add(fin.Get("TTB"))
        if "DY" in alistname: hist_template.Add(fin.Get("DY"))
        for i in range(1,hist_template.GetNbinsX()+1):
            if hist_template.GetBinContent(i) <=0 and i not in EMPTYX:
                EMPTYN+=1
                EMPTYX.append(i)
        EMPTYX.sort()
        if (EMPTYN)!=0: 
            print(EMPTYN, EMPTYX)
            print(binning)  
            for ib in range(EMPTYN):
                binning.pop(EMPTYX[ib*-1]-1)
            print('after pop: ', binning)
            #input()
        checkends=True
        if checkends:
            if len(binning)>2:
                if(binning[1]-binning[0])>2*(binning[2]-binning[1]): binning[0] = binning[1] - (binning[2]-binning[1])
            if len(binning)>4:
                if(binning[-1]-binning[-2])>2*(binning[-2]-binning[-3]): binning[-1] = binning[-2] + (binning[-2]-binning[-3])
            #print('modify the ends to make plot better', binning)
        fin.cd()
        # Copy all the hists
        for iprocs in alistname:
            hout[iprocs] = fin.Get(iprocs).Clone()
        fout = R.TFile.Open('{0}/{1}{2}'.format(args.out,addtional,ifile),'RECREATE')
        if 'mJ2' in ifile or EMPTYN==0:
            for iprocs in hout.keys():
                if checkends:
                    hsample[iprocs] = R.TH1F(iprocs,iprocs,len(binning)-1,array('d',binning))
                    for i in range(1,hout[iprocs].GetNbinsX()+1):
                            hsample[iprocs].SetBinContent(i,hout[iprocs].GetBinContent(i))
                            hsample[iprocs].SetBinError(i,hout[iprocs].GetBinError(i))
                    if 'CMS_vhh4b_scaleMuRMuF_ttH' in iprocs: hsample[iprocs].SetName(iprocs.replace('ttH_CMS',"ttH_hbb_CMS"))
                    elif 'ttH' in iprocs and 'ttH_hbb' not in iprocs:  hsample[iprocs].SetName(iprocs.replace("ttH",'ttH_hbb'))
                    hsample[iprocs].Write()
                else:
                    if 'CMS_vhh4b_scaleMuRMuF_ttH' in iprocs: hout[iprocs].SetName(iprocs.replace('ttH_CMS',"ttH_hbb_CMS"))
                    elif 'ttH' in iprocs and 'ttH_hbb' not in iprocs:  hout[iprocs].SetName(iprocs.replace("ttH",'ttH_hbb'))
                    hout[iprocs].Write()
            fout.Close()
        else:
            for iprocs in hout.keys():
                hsample[iprocs] = R.TH1F(iprocs,iprocs,len(binning)-1,array('d',binning))
                if iprocs in data_obs:
                    print("skip data")
                else:
                    ireal=0
                    for i in range(1,hout[iprocs].GetNbinsX()+1):
                        ireal+=1
                        if i in EMPTYX:
                            ireal-=1
                            merged = hsample[iprocs].GetBinContent(ireal) + hout[iprocs].GetBinContent(i)
                            merged_err = math.sqrt(hsample[iprocs].GetBinError(ireal)**2 + hout[iprocs].GetBinError(i)**2)
                            hsample[iprocs].SetBinContent(ireal,merged)
                            hsample[iprocs].SetBinError(ireal,merged_err)
                        else:
                            hsample[iprocs].SetBinContent(ireal,hout[iprocs].GetBinContent(i))
                            hsample[iprocs].SetBinError(ireal,hout[iprocs].GetBinError(i))
                if 'CMS_vhh4b_scaleMuRMuF_ttH' in iprocs: hsample[iprocs].SetName(iprocs.replace('ttH_CMS',"ttH_hbb_CMS"))
                elif 'ttH' in iprocs and 'ttH_hbb' not in iprocs:  hsample[iprocs].SetName(iprocs.replace("ttH",'ttH_hbb'))
                hsample[iprocs].Write()
            fout.Close()

def check_LHE():
    friend_LHE='test_221005_prepare_mergeLHE/'
    for ifile in VHHfiles:
        print(ifile)
        hsample = {}
        hout={}
        f1 = R.TFile.Open('{0}'.format(path+ifile),'READ')
        f2 = R.TFile.Open('{0}'.format(friend_LHE+ifile),'READ')
        alistname=['TT','TTB','ttH_hbb']
        for iprocs in alistname:
            print(iprocs)
            f1.cd()
            h1=f1.Get(iprocs).Clone()
            f2.cd()
            h2=f2.Get(iprocs).Clone()
            n1=h1.GetNbinsX()
            n2=h2.GetNbinsX()
            if n1!=n2:
                print('the number of bins are not equal! ', ifile)
                exit()
            else:
                for i in range(1,n1+1):
                    if h1.GetBinCenter(i)!=h2.GetBinCenter(i):
                        print('the Bin Contents are not equal! ', ifile, i, h1.GetBinContent(i), h2.GetBinContent(i))
                        exit()
                    if h1.GetBinContent(i)!=h2.GetBinContent(i):
                        print('the Bin Contents are not equal! ', ifile, i, h1.GetBinContent(i), h2.GetBinContent(i))
                        exit()
        print('finished, looks good')

def quickmodify():
    for ifile in VHHfiles:
        if 'Zll' not in ifile:
            continue
        hsample = {}
        hout={}
        fin = R.TFile.Open('{0}'.format(path+ifile),'READ')
        alist=fin.GetListOfKeys()
        alistname=[]
        for i in range(alist.GetEntries()):
            alistname.append(alist.At(i).GetName())
        fin.cd()
        for iprocs in alistname:
            hout[iprocs] = fin.Get(iprocs).Clone()

        fout = R.TFile.Open('{0}/{1}{2}'.format(args.out,addtional,ifile),'RECREATE')
        for iprocs in hout.keys():
            if '2L_R' in iprocs and 'TT' in iprocs: hout[iprocs].SetName(iprocs.replace('2L_R',"2L_tt"))
            elif '2L_R' in iprocs and 'DY' in iprocs: hout[iprocs].SetName(iprocs.replace('2L_R',"2L_dy"))
            hout[iprocs].Write()
        fout.Close()



#merge_signals()
#merge_years()
#merge_ZllTTCR_3b4b()
#merge_emptyBins()
quickmodify()
#
#check_LHE()




