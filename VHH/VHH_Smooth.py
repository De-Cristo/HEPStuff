# copied from https://gitlab.cern.ch/cms-hcg/ch-areas/VHLegacy/-/blob/master/scripts/modifyDC.py

import CombineHarvester.CombineTools.ch as ch
import ROOT
import time
import argparse
import re
import pprint
from collections import defaultdict
import numpy as np
from math import sqrt

def drop_zero_procs(chob,proc):
  print(chob,proc.name())
  null_yield = not (proc.rate() > 1e-6)
  if(null_yield):
    chob.FilterSysts(lambda sys: matching_proc(proc,sys))
  return null_yield


def drop_zero_systs(syst):
  null_yield = (not (syst.value_u() > 0. and syst.value_d()>0.) ) and syst.type() in 'shape'
  if(null_yield):
    print 'Dropping systematic ',syst.name(),' for region ', syst.bin(), ' ,process ', syst.process(), '. up norm is ', syst.value_u() , ' and down norm is ', syst.value_d()
  return null_yield


def remove_norm_effect(syst):
  print 'removing norm',syst
  syst.set_value_u(1.0)
  syst.set_value_d(1.0)

def symm(syst,nominal):
  print 'Symmetrising systematic ', syst.name(), ' in region ', syst.bin(), ' for process ', syst.process()
  hist_u = syst.ShapeUAsTH1F()
  hist_u.Scale(nominal.Integral()*syst.value_u())
  hist_d = nominal.Clone()
  hist_d.Scale(2)
  hist_d.Add(hist_u,-1)
  syst.set_shapes(hist_u,hist_d,nominal)

def smooth_markov(h):
    data = []
    data = [h.GetBinContent(i+1) for i in range(h.GetNBinsX())]
    s = ROOT.TSpectrum()
    smoothed = ROOT.TH1F("smoothed","smoothed",h.GetXaxis().GetXbins())
    s.SmoothMarkov(data,h.GetNBinsX(),3)
    for i in range(smoothed.GetNBinsX()): smoothed.SetBinContent(i+1,data[i])
    return smooth_markov

def smooth_lowess_tgraph(h):
    g=ROOT.TGraph()
    for i in range(h.GetNbinsX()): g.SetPoint(i,h.GetBinCenter(i+1), h.GetBinContent(i+1))
    smooth = ROOT.TGraphSmooth("normal")
    gout = smooth.SmoothLowess(g,"",0.5)
    hist_out = ROOT.TH1F('hist_out','hist_out',h.GetNbinsX(),h.GetXaxis().GetXbins().GetArray())
    for i in range(h.GetNbinsX()):
        x,y = ROOT.Double(0), ROOT.Double(0)
        gout.GetPoint(i,x,y)
        hist_out.Fill(x,y)
    return hist_out
    


def smoothshape(syst,proc):
        #print 'smoothing: ', syst, ' bin: ', syst.bin()
        nominal = proc.shape()
        nominal.Scale(proc.rate())
        # Up
        hist_u = syst.shape_u()
        hist_u.Scale(nominal.Integral()*syst.value_u())
        up_rate = nominal.Integral()*syst.value_u()
        hist_u.Divide(nominal)
        if nominal.GetNbinsX()>5:
          hist_u.Smooth(1)
        hist_u.Multiply(nominal)
        if hist_u.Integral()>0:hist_u.Scale(up_rate/hist_u.Integral())
        # Down
        hist_d = syst.shape_d()
        hist_d.Scale(nominal.Integral()*syst.value_d())
        dn_rate = nominal.Integral()*syst.value_d()
        hist_d.Divide(nominal)
        if nominal.GetNbinsX()>5:
          hist_d.Smooth(1)
        hist_d.Multiply(nominal)
        if hist_d.Integral()>0:hist_d.Scale(dn_rate/hist_d.Integral())
        # Set New shape
        if hist_u.Integral()>1e-5 and hist_d.Integral()>1e-5:
            syst.set_shapes(hist_u,hist_d,nominal)
        else:
            syst.set_shapes(nominal,nominal,nominal)

def matching_proc(p,s):
  return ((p.bin()==s.bin()) and (p.process()==s.process()) and (p.signal()==s.signal()) 
         and (p.analysis()==s.analysis()) and  (p.era()==s.era()) 
         and (p.channel()==s.channel()) and (p.bin_id()==s.bin_id()) and (p.mass()==s.mass()))

def symmetrise_syst(chob,proc,sys_name):
  nom_hist = proc.ShapeAsTH1F()
  nom_hist.Scale(proc.rate())
  chob.ForEachSyst(lambda s: symm(s,nom_hist) if (s.name()==sys_name and matching_proc(proc,s)) else None)

def symmetrise_smooth_syst(chob,syst):
  chob.cp().syst_name([syst.name()]).ForEachProc(lambda x: smoothshape(syst,x) if (matching_proc(x,syst)) else None)


def flipUpDown(syst,proc):
    #print 'flipUpDown: ', syst, ' bin: ', syst.bin()
    #print 'proc.shape(): ', proc.shape(), ' proc.rate() ', proc.rate()
    #print 'syst.value_u() ', syst.value_u(), ' syst.value_d() ', syst.value_d()
    nominal = proc.shape()
    hist_u = syst.shape_u()
    hist_d = syst.shape_d()
    name_u = hist_u.GetName()
    name_d = hist_d.GetName()
    #print name_u, name_d
    hist_u.SetTitle(name_d)
    hist_d.SetTitle(name_u)
    value_d = syst.value_d()
    value_u = syst.value_u()
    #print 'nominal.Integral() ', nominal.Integral(), ' hist_u.Integral() ', hist_u.Integral(), ' hist_d.Integral() ', hist_d.Integral()
    syst.set_shapes(hist_d,hist_u,nominal)
    syst.set_value_u(value_d)
    syst.set_value_d(value_u)

def flipUpDown_syst(chob,syst):
    chob.cp().syst_name([syst.name()]).ForEachProc(lambda x: flipUpDown(syst,x) if (matching_proc(x,syst)) else None)


def Checkshape(syst,proc):
    print 'proc.shape(): ', proc.shape(), ' proc.rate() ', proc.rate()
    print 'syst.value_u() ', syst.value_u(), ' syst.value_d() ', syst.value_d()

def Checkshape_syst(chob,syst):
    print('---------> ',syst.name())
    if( syst.value_u()<1e-4 or syst.value_d()<1e-4):
        print('syst.value_u() ', syst.value_u(), ' syst.value_d() ', syst.value_d(), syst.process())
        #chob.cp().FilterSysts(syst)
    #chob.cp().syst_name([syst.name()]).ForEachProc(lambda x: Checkshape(syst,x) if (matching_proc(x,syst)) else None)


def renameRatePara(chob,syst):
    #chob.cp().syst_name([syst.name()]).ForEachProc(lambda x: syst.set_name('test') if (matching_proc(x,syst) and syst.type()=='rateParam' and syst.process()=='TTB') else None)
    if syst.type()=='rateParam' and syst.process()=='TTB':
        name = syst.name()
        print name
        chob.cp().FilterSysts(lambda x: True if(x.type()=='rateParam' and x.process()=='TTB') else None)

    #    syst.set_name(name.replace('TT','TTB'))


def redefineBin(hist_data):
    oldbins = hist_data.GetXaxis().GetXbins().GetArray()
    oldbins_arr = np.array(oldbins)
    print(oldbins_arr)
    return oldbins_arr

parser = argparse.ArgumentParser()
parser.add_argument( '--channel', help='Comma separated list of channels to process (default=all)', default=None)
parser.add_argument( '--veto', help='Comma separated list of channels to veto (default=none)', default=None)
parser.add_argument( '--output_folder', help='Postfix string for output', default='')
parser.add_argument( '--input_folder', help='Directory with datacards', default='')
parser.add_argument( '--year', help='year', default='2017')
parser.add_argument( '--comb', help='running on comb datacards', default=False,action="store_true")

args = parser.parse_args()

output_folder = args.output_folder
year = args.year 
comb = args.comb 
select_chns = []
veto_chns = []
if args.channel is not None:
    select_chns = args.channel.split(',')
if args.veto is not None:
    veto_chns = args.veto.split(',')

ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')


input_folder = args.input_folder
if input_folder == '':
    print('please provide input_folder')
    exit()

all_cards = {
    'MET' : [
            input_folder+'/vhh4b_SB_Znn_101_13TeV'+year+'.txt',
            input_folder+'/vhh4b_SB_Znn_102_13TeV'+year+'.txt',
            input_folder+'/vhh4b_Znn_1_13TeV'+year+'.txt',
            input_folder+'/vhh4b_Znn_2_13TeV'+year+'.txt',
        ],
    'MET_Boosted' : [
            input_folder+'/vhh4b_SB_Znn_103_13TeV'+year+'.txt',
            input_folder+'/vhh4b_SB_Znn_104_13TeV'+year+'.txt',
            input_folder+'/vhh4b_Znn_3_13TeV'+year+'.txt',
            input_folder+'/vhh4b_Znn_4_13TeV'+year+'.txt',
        ],
    'SL' : [
            input_folder+'/vhh4b_Wln_5_13TeV'+year+'.txt',
            input_folder+'/vhh4b_Wln_6_13TeV'+year+'.txt',
            input_folder+'/vhh4b_SB_Wln_105_13TeV'+year+'.txt',
            input_folder+'/vhh4b_SB_Wln_106_13TeV'+year+'.txt',
        ],
    'SL_Boosted' : [
            input_folder+'/vhh4b_Wln_7_13TeV'+year+'.txt',
            input_folder+'/vhh4b_Wln_8_13TeV'+year+'.txt',
            input_folder+'/vhh4b_SB_Wln_107_13TeV'+year+'.txt',
            input_folder+'/vhh4b_SB_Wln_108_13TeV'+year+'.txt',
        ],

    'DL' : [
            input_folder+'/vhh4b_Zll_9_13TeV_run2.txt',
            input_folder+'/vhh4b_Zll_10_13TeV_run2.txt',
            input_folder+'/vhh4b_Zll_11_13TeV_run2.txt',
            input_folder+'/vhh4b_Zll_12_13TeV_run2.txt',
            input_folder+'/vhh4b_Zll_13_13TeV_run2.txt',
            input_folder+'/vhh4b_Zll_14_13TeV_run2.txt',
            input_folder+'/vhh4b_Zll_15_13TeV_run2.txt',
            input_folder+'/vhh4b_Zll_15_13TeV_run2.txt',
            input_folder+'/vhh4b_Zll_16_13TeV_run2.txt',
            input_folder+'/vhh4b_SB_Zll_109_13TeV_run2.txt',
            input_folder+'/vhh4b_SB_Zll_110_13TeV_run2.txt',
            input_folder+'/vhh4b_SB_Zll_111_13TeV_run2.txt',
        ],

}



if not comb:
    for chn in all_cards.keys():
        if (len(select_chns) > 0 and chn not in select_chns) or chn in veto_chns:
            print '>> Skipping %s' % chn
            continue
        print '>> Parsing %s cards' % chn
        for card in all_cards[chn]:
            cb = ch.CombineHarvester()
            cb.SetVerbosity(0)
            cb.ParseDatacard(card, analysis='vhh4b')
            #if str(year) == '2016':
                # Only flip for 2016
                #cb.cp().ForEachSyst(lambda x: flipUpDown_syst(cb,x) if (x.name()=='CMS_btag_LF_2016_2017_2018') else None)
            cb.cp().ForEachSyst(lambda x: symmetrise_smooth_syst(cb,x) if ( "CMS_scale_j" in x.name() or 'res_j' in x.name() or 'CMS_unclusteredEnergy' in x.name() or 'Reg' in x.name() or 'Pile' in x.name()) else None)
            cb.WriteDatacard(output_folder+'/'+card.split('/')[-1], output_folder+'/vhh4b_input_'+card.split('/')[-1].replace('.txt','.root'))


