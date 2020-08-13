import sys
import ROOT
#sys.path.append("/home/taohuang/DiHiggsAnalysis/CMSSW_9_4_0_pre1/src/HhhAnalysis/python/NanoAOD")
#sys.path.append("/home/taohuang/HhhAnalysis/python/NanoAOD")
import Samplelist as Slist
import os



faileddatasets =  ["WWToLNuQQ_aTGC_13TeV-madgraph-pythia8", "ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1"]
faileddatasets.append("ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1")

def get_event_weight_sum_file(filepath):
    tfile = ROOT.TFile(filepath, "READ")
    hist = tfile.Get("h_cutflow")
    event_weight_sum = hist.GetBinContent(1)
    tfile.Close()
    return event_weight_sum

def copy_hist_file1Tofile2(file1, file2, histname):
    #print "file1 ",file1, " file2 ",file2, " histname ",histname
    tfile1 = ROOT.TFile(file1, "READ")
    tfile2 = ROOT.TFile(file2, "UPDATE")
    histtmp = tfile2.Get(histname)
    if not histtmp:
        hist = tfile1.Get(histname)
        #print "hist ",hist
        histtmp = hist.Clone()
        hist.SetDirectory(tfile2)
        histtmp.Write()
    tfile2.Close()
    tfile1.Close()

def add_xsection_event_weight_sum(infile, xsection, event_weight_sum):
    tfile = ROOT.TFile(infile,"UPDATE")
    p = ROOT.TParameter(float)("cross_section", xsection)
    p2 = ROOT.TParameter(float)("event_weight_sum", event_weight_sum)
    p2.Write()
    p.Write()
    tfile.Close()

    #p = tfile.Get("cross_section")
    #p.SetVal( xsection )
def write_weight_MC(shortnames, dirWithLHEWeightScale, localSamplelist):
    for shortname in shortnames:
        for samplename in localSamplelist[shortname].keys():
            histname = "CountWeightedLHEWeightScale"
            histname2 = "h_cutflow"
            oldfile = os.path.join(dirWithLHEWeightScale, samplename+"_Friend.root")
            newfile = localSamplelist[shortname][samplename]["path"] 
            copy_hist_file1Tofile2(oldfile, newfile, histname)
            copy_hist_file1Tofile2(oldfile, newfile, histname2)
            add_xsection_event_weight_sum(newfile, localSamplelist[shortname][samplename]["cross_section"], localSamplelist[shortname][samplename]["event_weight_sum"])

def write_weight_MC_DY(shortnames, dirWithLHEWeightScale, localSamplelist, untagged_suffix):
    for shortname in shortnames:
        for samplename in localSamplelist[shortname+untagged_suffix].keys():
            histname = "CountWeightedLHEWeightScale"
            histname2 = "h_cutflow"
            oldsamplename = samplename.replace(untagged_suffix, "")
            oldfile = os.path.join(dirWithLHEWeightScale, oldsamplename+"_Friend.root")
            newfile = localSamplelist[shortname+untagged_suffix][samplename]["path"] 
            print "oldfile ",oldfile, " newfile ",newfile
            copy_hist_file1Tofile2(oldfile, newfile, histname)
            copy_hist_file1Tofile2(oldfile, newfile, histname2)
            add_xsection_event_weight_sum(newfile, localSamplelist[shortname+untagged_suffix][samplename]["cross_section"], localSamplelist[shortname+untagged_suffix][samplename]["event_weight_sum"])

            
            

    
full_local_samplelist = {}
#localdir = "/fdata/hepx/store/user/taohuang/HHNtuple_20180502_dataTT/"
datapath = "/Users/taohuang/Documents/DiHiggs/20180316_NanoAOD/"
localdir = datapath+"HHNtuple_20180518_addSys/"
localdir = datapath+"20180619_HHbbWW_addHME_addSys_10k/"
localdir = datapath+"20180619_HHbbWW_addHME_addSys_10k_addNN/"
addHME = True
addNN = True
for i,isample in enumerate(Slist.NumSample):
    
    if int(isample) <0:
        continue
    #if dataname in faileddatasets:
    #    continue
    sampleName = Slist.sampleN_short[i]
    #print "isample ",isample, " i ",i," Samplename ",sampleName, " dataset ",Slist.Nanodatasets[i]
    dataname =  Slist.Nanodatasets[i].split('/')[1]
    localfilepath = os.path.join(localdir, dataname+"_Friend.root")
    if addHME and addNN:
        localfilepath = os.path.join(localdir, dataname+"_HME_Friends_NN.root")
    if addHME and not(addNN):
        localfilepath = os.path.join(localdir, dataname+"_HME_Friends.root")
    if int(isample) < 13:
        localfilepath = os.path.join(datapath+"20180703_HHbbWW_addHME_addSys_10k_Signalonly/2018-07-05/", dataname+"_HME_Friends_NN.root")
    xsec =  Slist.MCxsections[i]
    if sampleName not in full_local_samplelist.keys():
        full_local_samplelist[sampleName] = {}

    full_local_samplelist[sampleName][dataname] = {}
    full_local_samplelist[sampleName][dataname]["path"] = localfilepath
    full_local_samplelist[sampleName][dataname]["cross_section"] = xsec
    #if sampleName == "TT":
        #full_local_samplelist[sampleName][dataname]["path"] = "HHNtuples_Run2017_20190306/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8_Friend.root"
    #    full_local_samplelist[sampleName][dataname]["path"] ="HHNtuples_Run2017_20190306/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_Run2017v2.root" 
    if not os.path.isfile(full_local_samplelist[sampleName][dataname]["path"]):
        print "error!! failed to find this file ", full_local_samplelist[sampleName][dataname]["path"]
    #print "to get event weight  sum ",localfilepath
    #full_local_samplelist[sampleName][dataname]["event_weight_sum"] =  get_event_weight_sum_file( localfilepath )
    

#print "full_local_samplelist ",full_local_samplelist
MCnames = ["TT","DY","sT","Wjet","VV", "ttV"]
signal = []
for mass in [260, 270, 300, 350, 400, 450, 500, 550, 600, 650,750, 800, 900]:
    signal.append('RadionM%d'%mass)
#write_weight_MC(MCnames, "HHNtuple_20180518_addSys/", full_local_samplelist)
#write_weight_MC(signal, "HHNtuple_20180518_addSys/", full_local_samplelist)

full_local_samplelist["Data"] = {}
datanames = ["DoubleMuon", "DoubleEG","MuonEG"]
for dataname in datanames:
    #localdatadir = "20180619_HHbbWW_addHME_addSys_10k/"
    localdatadir = localdir 
    full_local_samplelist["Data"][dataname] = {}
    full_local_samplelist["Data"][dataname]["path"] =  os.path.join(localdatadir,dataname+"Run2016.root")
    if addHME and addNN:
        full_local_samplelist["Data"][dataname]["path"] =  os.path.join(localdatadir,dataname+"_HME_Friends_NN.root")
    if addHME and not(addNN):
        full_local_samplelist["Data"][dataname]["path"] =  os.path.join(localdatadir,dataname+"_HME_Friends.root")


    if not os.path.isfile(full_local_samplelist["Data"][dataname]["path"]):
        print "error!! failed to find this file ", full_local_samplelist["Data"][dataname]["path"]


#print full_local_samplelist
untagged_MCname = ["TT", "sT"]
untagged_localdir = '20180619_HHbbWW_addHME_addSys_10k_DYestimation/'
untagged_localdir = localdir 
untagged_suffix = "_untagged"
for mcname in untagged_MCname:
    full_local_samplelist[mcname + untagged_suffix] = {}
    #for key in full_local_samplelist[mcname].keys():
    for key in full_local_samplelist[mcname].keys():
        dataname = key + untagged_suffix
        full_local_samplelist[mcname + untagged_suffix][key + untagged_suffix] = {}
        #full_local_samplelist[mcname][key]['path'] = os.path.join(untagged_localdir, key)
        full_local_samplelist[mcname + untagged_suffix][key + untagged_suffix]['path'] = os.path.join(untagged_localdir, key+"_Friend_untagged.root")
        if addHME and addNN:
            full_local_samplelist[mcname + untagged_suffix][key + untagged_suffix]['path'] = os.path.join(untagged_localdir, key+"_HME_addMbtagWeight_NN.root")
        if addHME and not(addNN):
            full_local_samplelist[mcname + untagged_suffix][key + untagged_suffix]['path'] = os.path.join(untagged_localdir, key+"_HME_addMbtagWeight.root")
        full_local_samplelist[mcname + untagged_suffix][key + untagged_suffix]['cross_section'] = full_local_samplelist[ mcname][key]['cross_section']
    if not os.path.isfile(full_local_samplelist[mcname + untagged_suffix][dataname]["path"]):
        print "error!! failed to find this file ", full_local_samplelist[mcname + untagged_suffix][dataname]["path"]
        #full_local_samplelist[mcname + untagged_suffix][key + untagged_suffix]['event_weight_sum'] =  get_event_weight_sum_file( full_local_samplelist[ mcname][key]['path'])

#write_weight_MC_DY(untagged_MCname, "HHNtuple_20180518_addSys/", full_local_samplelist, untagged_suffix)


full_local_samplelist["Data" + untagged_suffix] = {}
datas = ["DoubleMuon", "DoubleEG"]
untagged_localdatadir = untagged_localdir
for dataname in datas:
    full_local_samplelist["Data" + untagged_suffix][dataname] = {}
    full_local_samplelist["Data" + untagged_suffix][dataname]["path"] =  os.path.join(untagged_localdatadir,dataname+"_Run2016_untagged.root")
    if addHME and addNN:
        full_local_samplelist["Data" + untagged_suffix][dataname]["path"] =  os.path.join(untagged_localdatadir,dataname+"Run2016_HME_addMbtagWeight_NN.root")
    if addHME and not(addNN):
        full_local_samplelist["Data" + untagged_suffix][dataname]["path"] =  os.path.join(untagged_localdatadir,dataname+"Run2016_HME_addMbtagWeight.root")
    if not os.path.isfile(full_local_samplelist["Data" + untagged_suffix][dataname]["path"]):
        print "error!! failed to find this file ", full_local_samplelist["Data" + untagged_suffix][dataname]["path"]
#print "full_local_samplelist ",full_local_samplelist
