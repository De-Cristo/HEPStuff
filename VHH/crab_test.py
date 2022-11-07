import sys
from CRABClient.UserUtilities import config
config = config()

# V0X nanov7; V2X nanov8
thisversion = 'VHH4bPostNano-V22-2018'

#config.General.requestName = 'VHH4bPostNano2018_MC_V13_p02'
config.General.workArea = '/afs/cern.ch/work/l/lichengz/VHH/crabarea/'

config.General.transferOutputs = True
config.General.transferLogs = False # cern cms eos is limited by 10K files

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script.sh'
config.JobType.scriptArgs = ['isMC=1','era=UL2018','dataRun=X','isVjets=0']
config.JobType.inputFiles = ['../keep_and_drop.txt','../postproc.py','../../../../../../scripts/haddnano.py','../Spring16_25nsV6_DATA_PtResolution_AK4PFchs.txt'] #hadd nano will not be needed once nano tools are in cmssw
config.JobType.sendPythonFolder	 = True

#config.Data.inputDataset = '/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2017_TrancheIV_v6-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 200000 # too large => exceed wall-time/disk; too small => file too many, exceed eos file number
#config.Data.totalUnits = 10
config.Data.outLFNDirBase = '/store/group/phys_higgs/hbb/ntuples/VHH4bPostNanov8/{0}/'.format(thisversion)
config.Data.publication = True
#config.Data.outputDatasetTag = 'RunIISummer16MiniAODv2-PUMoriond17-80X-VHbbPostNano2017_V1'
config.Data.allowNonValidInputDataset = True
#config.Site.storageSite = 'T3_CH_CERNBOX'
#config.Site.storageSite = 'T2_CN_Beijing'
config.Site.storageSite = 'T2_CH_CERN'

config.JobType.maxMemoryMB = 4000

#sites=['T2_CH_CERN', 'T2_DE_DESY', 'T2_US_FNAL']

if __name__ == '__main__':
    f=open(sys.argv[1]) 
    content = f.readlines()
    content = [x.strip() for x in content] 
    from CRABAPI.RawCommand import crabCommand
    n=700
    for dataset in content :
	#site=sites[n%4]
	#config.Site.storageSite=site
	#if site=='T2_CH_CERN' :
	#	config.Data.outLFNDirBase=  '/store/group/cmst3/group/nanoAOD/NanoTestProd006'
	#else :
        #	config.Data.outLFNDirBase = '/store/user/%s/NanoTestProd006/' % (getUsernameFromSiteDB())

        config.Data.inputDataset = dataset
	n+=1
	nnn="%s"%n
        #config.General.requestName = thisversion+"_"+dataset.split('/')[1].split('-')[0]+dataset.split('/')[2].split('-')[0]+'_'+nnn # name gets too long [forbiden by crab]
        config.General.requestName = '_'.join(thisversion.split('-')[1:])+"_"+dataset.split('/')[1].split('-')[0]+dataset.split('/')[2].split('-')[0]+'_'+nnn
        config.Data.outputDatasetTag = thisversion+dataset.split('/')[2].split('-')[0]
        crabCommand('submit', config = config)

