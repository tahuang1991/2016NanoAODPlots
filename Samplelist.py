#!/usr/bin/python
import os
import sys

fdatadir = "/fdata/hepx/store/user/taohuang/NANOAOD/"
datasets  = []; NumSample = []; sampleN_short = []
Nanodatasets = []; localdirs = {}
MCxsections = []
#doTT=True; doDY=True; doVV=True; doSingleT=True; doWjets=True; dottV=True

##DoubleEG
datasets.append('/DoubleEG/Run2016B-05Feb2018_ver1-v1/NANOAOD')
NumSample.append('-1'); sampleN_short.append('DoubleEGRun2016Bver1')
MCxsections.append(-1.0)
datasets.append('/DoubleEG/Run2016B-05Feb2018_ver2-v1/NANOAOD')
NumSample.append('-2'); sampleN_short.append('DoubleEGRun2016Bver2')
MCxsections.append(-1.0)
datasets.append('/DoubleEG/Run2016C-05Feb2018-v1/NANOAOD')
NumSample.append('-3'); sampleN_short.append('DoubleEGRun2016C')
MCxsections.append(-1.0)
datasets.append('/DoubleEG/Run2016D-05Feb2018-v1/NANOAOD')
NumSample.append('-4'); sampleN_short.append('DoubleEGRun2016D')
MCxsections.append(-1.0)
datasets.append('/DoubleEG/Run2016E-05Feb2018-v1/NANOAOD')
NumSample.append('-5'); sampleN_short.append('DoubleEGRun2016E')
MCxsections.append(-1.0)
datasets.append('/DoubleEG/Run2016F-05Feb2018-v1/NANOAOD')
NumSample.append('-6'); sampleN_short.append('DoubleEGRun2016F')
MCxsections.append(-1.0)
datasets.append('/DoubleEG/Run2016G-05Feb2018-v1/NANOAOD')
NumSample.append('-7'); sampleN_short.append('DoubleEGRun2016G')
MCxsections.append(-1.0)
datasets.append('/DoubleEG/Run2016H-05Feb2018_ver2-v1/NANOAOD')
NumSample.append('-8'); sampleN_short.append('DoubleEGRun2016Hver2')
MCxsections.append(-1.0)
datasets.append('/DoubleEG/Run2016H-05Feb2018_ver3-v1/NANOAOD')
NumSample.append('-9'); sampleN_short.append('DoubleEGRun2016Hver3')
MCxsections.append(-1.0)
##DoubleMuon
datasets.append('/DoubleMuon/Run2016B-05Feb2018_ver1-v1/NANOAOD')
NumSample.append('-10'); sampleN_short.append('DoubleMuonRun2016Bver1')
MCxsections.append(-1.0)
datasets.append('/DoubleMuon/Run2016B-05Feb2018_ver2-v1/NANOAOD')
NumSample.append('-11'); sampleN_short.append('DoubleMuonRun2016Bver2')
MCxsections.append(-1.0)
datasets.append('/DoubleMuon/Run2016C-05Feb2018-v1/NANOAOD')
NumSample.append('-12'); sampleN_short.append('DoubleMuonRun2016C')
MCxsections.append(-1.0)
datasets.append('/DoubleMuon/Run2016D-05Feb2018-v1/NANOAOD')
NumSample.append('-13'); sampleN_short.append('DoubleMuonRun2016D')
MCxsections.append(-1.0)
datasets.append('/DoubleMuon/Run2016E-05Feb2018-v1/NANOAOD')
NumSample.append('-14'); sampleN_short.append('DoubleMuonRun2016E')
MCxsections.append(-1.0)
datasets.append('/DoubleMuon/Run2016F-05Feb2018-v1/NANOAOD')
NumSample.append('-15'); sampleN_short.append('DoubleMuonRun2016F')
MCxsections.append(-1.0)
datasets.append('/DoubleMuon/Run2016G-05Feb2018-v1/NANOAOD')
NumSample.append('-16'); sampleN_short.append('DoubleMuonRun2016G')
MCxsections.append(-1.0)
datasets.append('/DoubleMuon/Run2016H-05Feb2018_ver2-v1/NANOAOD')
NumSample.append('-17'); sampleN_short.append('DoubleMuonRun2016Hver2')
MCxsections.append(-1.0)
datasets.append('/DoubleMuon/Run2016H-05Feb2018_ver3-v1/NANOAOD')
NumSample.append('-18'); sampleN_short.append('DoubleMuonRun2016Hver3')
MCxsections.append(-1.0)
#MuonEG
datasets.append('/MuonEG/Run2016B-05Feb2018_ver1-v1/NANOAOD')
NumSample.append('-19'); sampleN_short.append('MuonEGRun2016Bver2')
MCxsections.append(-1.0)
datasets.append('/MuonEG/Run2016B-05Feb2018_ver2-v1/NANOAOD')
NumSample.append('-20'); sampleN_short.append('MuonEGRun2016Bver2')
MCxsections.append(-1.0)
datasets.append('/MuonEG/Run2016C-05Feb2018-v1/NANOAOD')
NumSample.append('-21'); sampleN_short.append('MuonEGRun2016C')
MCxsections.append(-1.0)
datasets.append('/MuonEG/Run2016D-05Feb2018-v1/NANOAOD')
NumSample.append('-22'); sampleN_short.append('MuonEGRun2016D')
MCxsections.append(-1.0)
datasets.append('/MuonEG/Run2016E-05Feb2018-v1/NANOAOD')
NumSample.append('-23'); sampleN_short.append('MuonEGRun2016E')
MCxsections.append(-1.0)
datasets.append('/MuonEG/Run2016F-05Feb2018-v1/NANOAOD')
NumSample.append('-24'); sampleN_short.append('MuonEGRun2016F')
MCxsections.append(-1.0)
datasets.append('/MuonEG/Run2016G-05Feb2018-v1/NANOAOD')
NumSample.append('-25'); sampleN_short.append('MuonEGRun2016G')
MCxsections.append(-1.0)
datasets.append('/MuonEG/Run2016H-05Feb2018_ver2-v1/NANOAOD')
NumSample.append('-26'); sampleN_short.append('MuonEGRun2016Hver2')
MCxsections.append(-1.0)
datasets.append('/MuonEG/Run2016H-05Feb2018_ver3-v1/NANOAOD')
NumSample.append('-27'); sampleN_short.append('MuonEGRun2016Hver3')
MCxsections.append(-1.0)


masspoints = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900]
for mass in masspoints:
    datasets.append("/GluGluToRadionToHHTo2B2VTo2L2Nu_M-%d_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM"%mass)
    NumSample.append(masspoints.index(mass)); sampleN_short.append('RadionM%d'%mass)
    MCxsections.append(5.0)#by default, assume the cross section for signal is 5pb
#datasets.append("/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-*_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
#NumSample.append('2'); sampleN_short.append('Graviton')

# TT## FIXME, use official one later
#datasets.append('/TTTo2L2Nu_13TeV-powheg/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/MINIAODSIM')
datasets.append('/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
#datasets.append('/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/arizzi-RunIIFall17MiniAOD-94X-Nano01Fall17-e273b12d9f89d622a34e4bc98b05ee29/USER')
NumSample.append('13'); sampleN_short.append('TT')
#MCxsections.append(72.1)
#MCxsections.append(76.7)
MCxsections.append(87.31)
# DY
#datasets.append('/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')

datasets.append('/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('14'); sampleN_short.append('DY')
MCxsections.append(18610.0)
datasets.append('/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
NumSample.append('15'); sampleN_short.append('DY')
MCxsections.append(4758.9)
datasets.append('/DYToLL_1J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
NumSample.append('16'); sampleN_short.append('DY')
MCxsections.append(929.1)
datasets.append('/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
NumSample.append('17'); sampleN_short.append('DY')
MCxsections.append(337.1)
# VV
datasets.append('/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('18'); sampleN_short.append('VV')
MCxsections.append(3.22)
datasets.append('/ZZTo2L2Nu_13TeV_powheg_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('19'); sampleN_short.append('VV')
MCxsections.append(0.564)
datasets.append('/ZZTo4L_13TeV_powheg_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('20'); sampleN_short.append('VV')
MCxsections.append(1.256)
#datasets.append('/WWToLNuQQ_aTGC_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
#NumSample.append('21'); sampleN_short.append('VV')
#MCxsections.append(49.997)# ## not available now because of pdf uncertainty
#FIXME
#datasets.append('/WWTo2L2Nu_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
#datasets.append('/WWTo2L2Nu_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
#NumSample.append('22'); sampleN_short.append('VV') ### not available now
#MCxsections.append(12.178)
datasets.append('/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('23'); sampleN_short.append('VV')
MCxsections.append(5.595)
#FIXME
#datasets.append('/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
#NumSample.append('24'); sampleN_short.append('VV') ### not available now 
#MCxsections.append(3.033)
datasets.append('/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3/MINIAODSIM')
NumSample.append('25'); sampleN_short.append('VV')
MCxsections.append(10.71)
datasets.append('/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('26'); sampleN_short.append('VV')
MCxsections.append(4.42965)
##sT
datasets.append('/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('27'); sampleN_short.append('sT')
MCxsections.append(136.02)
datasets.append('/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('28'); sampleN_short.append('sT')
MCxsections.append(80.95)
datasets.append('/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('29'); sampleN_short.append('sT')
MCxsections.append(3.36)
datasets.append('/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
NumSample.append('30'); sampleN_short.append('sT')
MCxsections.append(19.5545)
datasets.append('/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
NumSample.append('31'); sampleN_short.append('sT')
MCxsections.append(19.5545)
# W + Jets
datasets.append('/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('32'); sampleN_short.append('Wjet')
MCxsections.append(61526.7)
datasets.append('/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM')
NumSample.append('33'); sampleN_short.append('Wjet')
MCxsections.append(1627.45)

datasets.append('/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM')
NumSample.append('34'); sampleN_short.append('Wjet')
MCxsections.append(435.237)
datasets.append('/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
NumSample.append('35'); sampleN_short.append('Wjet')
MCxsections.append(59.181)
#FIXME
#datasets.append('/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
#NumSample.append('36'); sampleN_short.append('Wjet')### not available now
MCxsections.append(14.58)
datasets.append('/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
NumSample.append('37'); sampleN_short.append('Wjet')
MCxsections.append(6.656)
datasets.append('/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
NumSample.append('38'); sampleN_short.append('Wjet')
MCxsections.append(1.608)
datasets.append('/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
NumSample.append('39'); sampleN_short.append('Wjet')
MCxsections.append(0.0389)
# tt + V
datasets.append('/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('40'); sampleN_short.append('ttV')
MCxsections.append(0.4062)
datasets.append('/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM')
NumSample.append('41'); sampleN_short.append('ttV')
MCxsections.append(0.2043)
datasets.append('/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('42'); sampleN_short.append('ttV')
MCxsections.append(0.5297)
datasets.append('/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext3-v1/MINIAODSIM')
NumSample.append('43'); sampleN_short.append('ttV')
MCxsections.append(0.2529)

alljobtypes = set(sampleN_short)
for job in alljobtypes:
    localdirs[job] = []

for ijob, job in enumerate(datasets):
    nsample = int(NumSample[ijob])
    jobtype = sampleN_short[ijob]
    dataname = ""
    datadir = " "
    #print "nsample ",nsample, " jobtype ",jobtype
    if nsample < 0:
        datadir = sampleN_short[ijob]
	dataname = job
        #print "real data nsample ",nsample, " datadir ",datadir
    elif nsample > 0:
        datadir = job.split('/')[1]
        #print "MC nsample ",nsample, " datadir ",datadir, "MiniAOD dataset ",job.split('/')
	#query = "dataset dataset=/%s/*/NANOAODSIM"%(datadir)
        #pdata = os.popen("dasgoclient -limit=0 -query='{query}'".format(query = query))	
        #founddataset = False
	#for line in pdata:
	#    #print "dataset ",line," datatype ",datadir
	#    if datadir in line:
	#        founddataset = True
	#        dataname = line[:-1]	
	#if not(founddataset): 
	#    print "WARNING!!!!! no dataset found for ",datadir
    localdirs[jobtype].append(os.path.join(fdatadir, datadir))

Nanodatasets.append("/DoubleEG/Run2016B-05Feb2018_ver1-v1/NANOAOD")
Nanodatasets.append("/DoubleEG/Run2016B-05Feb2018_ver2-v1/NANOAOD")
Nanodatasets.append("/DoubleEG/Run2016C-05Feb2018-v1/NANOAOD")
Nanodatasets.append("/DoubleEG/Run2016D-05Feb2018-v1/NANOAOD")
Nanodatasets.append("/DoubleEG/Run2016E-05Feb2018-v1/NANOAOD")
Nanodatasets.append("/DoubleEG/Run2016F-05Feb2018-v1/NANOAOD")
Nanodatasets.append("/DoubleEG/Run2016G-05Feb2018-v1/NANOAOD")
Nanodatasets.append("/DoubleEG/Run2016H-05Feb2018_ver2-v1/NANOAOD")
Nanodatasets.append("/DoubleEG/Run2016H-05Feb2018_ver3-v1/NANOAOD")
Nanodatasets.append("/DoubleMuon/Run2016B-05Feb2018_ver1-v1/NANOAOD")
Nanodatasets.append("/DoubleMuon/Run2016B-05Feb2018_ver2-v1/NANOAOD")
Nanodatasets.append("/DoubleMuon/Run2016C-05Feb2018-v1/NANOAOD")
Nanodatasets.append("/DoubleMuon/Run2016D-05Feb2018-v1/NANOAOD")
Nanodatasets.append("/DoubleMuon/Run2016E-05Feb2018-v1/NANOAOD")
Nanodatasets.append("/DoubleMuon/Run2016F-05Feb2018-v1/NANOAOD")
Nanodatasets.append("/DoubleMuon/Run2016G-05Feb2018-v1/NANOAOD")
Nanodatasets.append("/DoubleMuon/Run2016H-05Feb2018_ver2-v1/NANOAOD")
Nanodatasets.append("/DoubleMuon/Run2016H-05Feb2018_ver3-v1/NANOAOD")
Nanodatasets.append("/MuonEG/Run2016B-05Feb2018_ver1-v1/NANOAOD")
Nanodatasets.append("/MuonEG/Run2016B-05Feb2018_ver2-v1/NANOAOD")
Nanodatasets.append("/MuonEG/Run2016C-05Feb2018-v1/NANOAOD")
Nanodatasets.append("/MuonEG/Run2016D-05Feb2018-v1/NANOAOD")
Nanodatasets.append("/MuonEG/Run2016E-05Feb2018-v1/NANOAOD")
Nanodatasets.append("/MuonEG/Run2016F-05Feb2018-v1/NANOAOD")
Nanodatasets.append("/MuonEG/Run2016G-05Feb2018-v1/NANOAOD")
Nanodatasets.append("/MuonEG/Run2016H-05Feb2018_ver2-v1/NANOAOD")
Nanodatasets.append("/MuonEG/Run2016H-05Feb2018_ver3-v1/NANOAOD")
Nanodatasets.append("/GluGluToRadionToHHTo2B2VTo2L2Nu_M-260_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/GluGluToRadionToHHTo2B2VTo2L2Nu_M-270_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/GluGluToRadionToHHTo2B2VTo2L2Nu_M-300_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/GluGluToRadionToHHTo2B2VTo2L2Nu_M-350_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/GluGluToRadionToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/GluGluToRadionToHHTo2B2VTo2L2Nu_M-450_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/GluGluToRadionToHHTo2B2VTo2L2Nu_M-500_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/GluGluToRadionToHHTo2B2VTo2L2Nu_M-550_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/GluGluToRadionToHHTo2B2VTo2L2Nu_M-600_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/GluGluToRadionToHHTo2B2VTo2L2Nu_M-650_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/GluGluToRadionToHHTo2B2VTo2L2Nu_M-750_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/GluGluToRadionToHHTo2B2VTo2L2Nu_M-800_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/GluGluToRadionToHHTo2B2VTo2L2Nu_M-900_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
#TTbar
#Nanodatasets.append("/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/arizzi-RunIIFall17MiniAOD-94X-Nano01Fall17-e273b12d9f89d622a34e4bc98b05ee29/USER")
Nanodatasets.append('/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM')

# DY
Nanodatasets.append("/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM")
Nanodatasets.append("/DYToLL_1J_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM")
Nanodatasets.append("/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM")
# VV
Nanodatasets.append("/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/ZZTo2L2Nu_13TeV_powheg_pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/ZZTo4L_13TeV_powheg_pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
#Nanodatasets.append("/WWToLNuQQ_aTGC_13TeV-madgraph-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM")
#sT
Nanodatasets.append("/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
#W+jets
Nanodatasets.append("/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext2-v1/NANOAODSIM")
Nanodatasets.append("/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext2-v1/NANOAODSIM")
Nanodatasets.append("/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext2-v1/NANOAODSIM")
Nanodatasets.append("/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM")
Nanodatasets.append("/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM")
Nanodatasets.append("/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM")
# tt + V
Nanodatasets.append("/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext2-v1/NANOAODSIM")
Nanodatasets.append("/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
Nanodatasets.append("/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext3-v1/NANOAODSIM")



outAnalist = {}
outAnadir = "/fdata/hepx/store/user/taohuang/HHNtuple_20180328_fixedleptonDZeff/"
for i,datasetname in enumerate( Nanodatasets ):
    sampleName = sampleN_short[i]
    if NumSample[i] < 0:
    	sampleName = "Data"
    outAnafile = os.path.join(outAnadir, Nanodatasets[i].split('/')[1])
    if hasattr(outAnalist, sampleName):
	outAnalist[sampleName].append(outAnafile)
    else:
	outAnalist[sampleName] = []
	outAnalist[sampleName].append(outAnafile)

dataintxt = open("2016MCSamplelist.txt","w+")
for dataset in datasets:
    dataintxt.write(dataset+"\n")
dataintxt.close()
