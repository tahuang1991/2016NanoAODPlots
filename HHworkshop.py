import ROOT
import os
from localSamplelist import * 
from bbWWPlotterSystematics import *
import numpy as np
from math import sqrt

import sys
sys.argv.append( '-b' )
sys.argv.append( '-q' )

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)


ROOT.gStyle.SetOptStat(0)

TotalLumi = 36.8 #fb-1
#TotalLumi = 35.920
#extraweight is used for compensating the "HLT" safe ID cut, for Louvain ntuples
channelcuts = {"MuMu":{"cut":"isMuMu","extraweight": 1.0,"Data":"DoubleMuon", "latex":"#mu#mu"},
        "MuEl":{"cut":"(isMuEl || isElMu)", "extraweight": 1.0,"Data":"MuonEG", "latex":"e#mu"},
        "ElEl":{"cut":"isElEl","extraweight":1.0,"Data":"DoubleEG", "latex":"ee"},
		}
#cutflows = ["All","Trigger","online-offline matching","dilepton PtEta","dilepton IP","dilepton ID","dilepton Iso","HLT Safe ID","nlepton>=3 veto","M_{ll}>12","NJets>=2","dijet PtEta","DR_j_l > 0.3","dijet btagging","M_{ll}<76"]
#cutflows = ["All","dilepton PtEta","dilepton IP","dilepton ID","dilepton Iso","HLT Safe ID","EMTFBug","HLT matching","M_{ll}>12","NJets>=2","dijet PtEta","DR_j_l > 0.3","dijet btagging","M_{ll}<76"]
cutflows = ["All","dilepton PtEta","dilepton IP","dilepton ID","dilepton Iso","HLT Safe ID","HLT matching","EMTFBUG","M_{ll}>12","NJets>=2","dijet PtEta","DR_j_l > 0.3","dijet btagging","M_{ll}<76"]
def get_xsection(shortname, samplename = ''):
    if len(full_local_samplelist[shortname].keys()) == 1:
        samplename = full_local_samplelist[shortname].keys()[0]
    elif len(full_local_samplelist[shortname].keys()) > 1 and samplename == '':
        raise ValueError("no proper samplename found: ",samplename )

    return full_local_samplelist[shortname][samplename]["cross_section"]
def get_xsection_file(filename):
    tfile = ROOT.TFile(filename, "READ")
    xsec = tfile.Get("cross_section")
    return xsec.GetVal()

def get_event_weight_sum_file(filepath):
    tfile = ROOT.TFile(filepath, "READ")
    hist = tfile.Get("h_cutflow")
    event_weight_sum = hist.GetBinContent(1)
    tfile.Close()
    return event_weight_sum

def get_event_weight_sum(shortname, samplename=''):
    if len(full_local_samplelist[shortname].keys()) == 1:
        samplename = full_local_samplelist[shortname].keys()[0]
    elif len(full_local_samplelist[shortname].keys()) > 1 and samplename == '':
        raise ValueError("no proper samplename found: ",samplename )

    filepath = full_local_samplelist[shortname][samplename]["path"]
    #print "samplename ",samplename, " file ",filepath
    return get_event_weight_sum_file(filepath)

def plotCutflowHist_data(outdir):
    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2]
    channels = ["ElEl","MuEl","MuMu"]
    #channels = ["ElEl","MuMu"]

    legend = ROOT.TLegend(0.65,0.62,0.88,0.84); 
    legend.SetTextSize(0.04); legend.SetTextFont(42)
    legend.SetHeader("cutflow")
    #legend.SetBorderSize(0)
    hs = ROOT.THStack("hs"," ")
    allhist = []
    shortname = "Data"
    for i,ch in enumerate(channels):
        dataname = channelcuts[ch]["Data"]
	filepath = full_local_samplelist[shortname][dataname]["path"]
	tfile = ROOT.TFile(filepath, "READ")
	print "samplename ",dataname, " file ",filepath," h_cutflow_"+dataname
        allhist.append(tfile.Get("h_cutflow_"+dataname))
        allhist[-1].SetDirectory(0) 
        allhist[-1].SetLineColor(colors[i])
        allhist[-1].SetLineWidth(2)
        ch_yield = allhist[-1].GetBinContent(len(cutflows))
        hs.Add(allhist[-1])
        entry = legend.AddEntry(allhist[-1], channelcuts[ch]["latex"]+": %.1f"%ch_yield,"l")
        entry.SetTextColor(colors[i])
        tfile.Close()

    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetLogy()
    tex0 = ROOT.TLatex(0.1,0.92, " #scale[1.4]{#font[61]{CMS}} Internal: Work in pogress"+"  "*30+"35.87 fb^{-1} (13 TeV),2016")
    tex0.SetNDC(); tex0.SetTextSize(.045); tex0.SetTextFont(42)
    tex1 = ROOT.TLatex(0.15, 0.83, "Data, Run2016 ")
    tex1.SetNDC(); tex1.SetTextSize(.035); tex1.SetTextFont(42)
    hs.Draw("nostackhist")
    xaxis = hs.GetHistogram().GetXaxis()
    for i, cut in enumerate(cutflows):
        xaxis.SetBinLabel(i+1, cut)
    legend.Draw("same")
    tex0.Draw("same")
    tex1.Draw("same")
    c1.SaveAs(outdir+"Data_Run2016_cutflow.pdf")

def plotCutflowHist(outdir, shortname, samplename = ""):
    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2]
    xsec = 1.0; event_weight_sum = 1.0
    if shortname != "Data":
	xsec = get_xsection(shortname, samplename)
	event_weight_sum = get_event_weight_sum(shortname, samplename)
    channels = ["ElEl","MuEl","MuMu"]
    if len(full_local_samplelist[shortname].keys()) == 1:
        samplename = full_local_samplelist[shortname].keys()[0]
    elif len(full_local_samplelist[shortname].keys()) > 1 and samplename == '':
        raise ValueError("no proper samplename found: ",samplename )

    filepath = full_local_samplelist[shortname][samplename]["path"]
    #print "samplename ",samplename, " file ",filepath
    tfile = ROOT.TFile(filepath, "READ")
    legend = ROOT.TLegend(0.65,0.62,0.88,0.84); 
    legend.SetTextSize(0.04); legend.SetTextFont(42)
    legend.SetHeader("cutflow")
    #legend.SetBorderSize(0)
    hs = ROOT.THStack("hs"," ")
    allhist = []
    for i,ch in enumerate(channels):
        datasetname = channelcuts[ch]["Data"]
        allhist.append(tfile.Get("h_cutflow_"+datasetname))
        weight = TotalLumi*xsec*1000.0/event_weight_sum
	if shortname == "Data":
	       weight = 1.0
        allhist[-1].Scale(weight)
        allhist[-1].SetLineColor(colors[i])
        allhist[-1].SetLineWidth(2)
        ch_yield = allhist[-1].GetBinContent(len(cutflows))
        hs.Add(allhist[-1])
        entry = legend.AddEntry(allhist[-1], channelcuts[ch]["latex"]+": %.1f"%ch_yield,"l")
        entry.SetTextColor(colors[i])

    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetLogy()
    tex0 = ROOT.TLatex(0.1,0.92, " #scale[1.4]{#font[61]{CMS}} Work in progress"+"  "*30+"35.87 fb^{-1} (13 TeV),2016")
    tex0.SetNDC(); tex0.SetTextSize(.045); tex0.SetTextFont(42)
    tex1 = ROOT.TLatex(0.15, 0.83, "MC: "+samplename)
    tex1.SetNDC(); tex1.SetTextSize(.035); tex1.SetTextFont(42)
    hs.Draw("nostackhist")
    xaxis = hs.GetHistogram().GetXaxis()
    for i, cut in enumerate(cutflows):
        xaxis.SetBinLabel(i+1, cut)
    legend.Draw("same")
    tex0.Draw("same")
    tex1.Draw("same")
    c1.SaveAs(outdir+samplename+"_Run2016_cutflow.pdf")


def runallCutflowhist(outdir, mcnames):
    for i, key in enumerate(mcnames):
        shortname = key 
        for iname, samplename in enumerate(full_local_samplelist[key].keys()):
            plotCutflowHist(outdir, shortname, samplename )


def plotCutflowHist_allMC(outdir, bgnames):
    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2]
    channels = ["ElEl","MuEl","MuMu"]

    #print "samplename ",samplename, " file ",filepath
    legend = ROOT.TLegend(0.62,0.62,0.88,0.84); 
    legend.SetTextSize(0.04); legend.SetTextFont(42)
    legend.SetHeader("all backgrouds cutflow")
    #legend.SetBorderSize(0)
    hs = ROOT.THStack("hs"," ")
    allhist = []
    nbins = len(cutflows)+4
    for ich,ch in enumerate(channels):
        datasetname = channelcuts[ch]["Data"]
        ch_hist = ROOT.TH1F("cutflow_"+ch, "cutflow_"+ch, nbins, 0, nbins)
	for i, key in enumerate(bgnames):
            shortname = key 
            for iname, samplename in enumerate(full_local_samplelist[key].keys()):
                filepath = full_local_samplelist[shortname][samplename]["path"]
                tfile = ROOT.TFile(filepath, "READ")
                xsec = get_xsection(shortname, samplename)
                event_weight_sum = get_event_weight_sum(shortname, samplename)
                hist = tfile.Get("h_cutflow_"+datasetname)
                weight = TotalLumi*xsec*1000.0/event_weight_sum
                hist.Scale(weight)
                ch_hist.Add(hist)
        allhist.append(ch_hist)
        allhist[-1].SetLineColor(colors[ich])
        allhist[-1].SetLineWidth(2)
        ch_yield = allhist[-1].GetBinContent(len(cutflows))
        hs.Add(allhist[-1])
        entry = legend.AddEntry(allhist[-1], channelcuts[ch]["latex"]+": %.1f"%ch_yield,"l")
        entry.SetTextColor(colors[ich])

    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetLogy()
    tex0 = ROOT.TLatex(0.1,0.92, " #scale[1.4]{#font[61]{CMS}} Work in progress"+"  "*30+"35.87 fb^{-1} (13 TeV),2016")
    tex0.SetNDC(); tex0.SetTextSize(.045); tex0.SetTextFont(42)
    hs.Draw("nostackhist")
    xaxis = hs.GetHistogram().GetXaxis()
    for i, cut in enumerate(cutflows):
        xaxis.SetBinLabel(i+1, cut)
    legend.Draw("same")
    tex0.Draw("same")
    c1.SaveAs(outdir+"HHbbWW_backgrounds_Run2016_cutflow.pdf")


def DrellYanDataDriven(channel, filedict, todraw, cut, xbins, xtitle, suffix, makeDYplots, plotname):
    treename = "Friends"

    if len(xbins) == 3:
        nbins = xbins[0]; xmin = xbins[1]; xmax =  xbins[2]
        xbins = []
        binwidth = (xmax-xmin)/nbins
        for i in range(0, nbins+1):
            xbins.append(xmin + i*binwidth)
        xbins = np.asarray(xbins)

    untagged_suffix = "_untagged"
    hist_data = ROOT.TH1F("data" + untagged_suffix +"_"+channel+"_%s"%(suffix), "Data" + untagged_suffix+"_"+channel+"_%s"%(suffix), len(xbins)-1, xbins)
    hist_TT = ROOT.TH1F("TT"+ untagged_suffix+"_"+channel+"_%s"%(suffix), "TT"+ untagged_suffix+"_"+channel+"_%s"%(suffix), len(xbins)-1, xbins)
    TT_dict = filedict["TT" + untagged_suffix]["TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8"+ untagged_suffix]
    fileTT = TT_dict['path']
    xsec = TT_dict['cross_section']
    #xsec = get_xsection_file(fileTT)
    event_weight_sum = get_event_weight_sum_file(full_local_samplelist['TT']['TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8']['path'])
    Mbtag_weight = "dy_Mbtag_weight"
    #Mbtag_weight = "1"
    weight = Mbtag_weight+"*sample_weight*event_reco_weight*{totallumi}*{cross_section}*1000.0/{event_weight_sum}".format(totallumi = TotalLumi, cross_section = xsec, event_weight_sum = event_weight_sum)
    #print "channel ",channel
    #print "TT cross_section ",xsec, " event_weight_sum ",event_weight_sum
    finalcut = "("+ cut + " && "+ channelcuts[channel]["cut"] +")*"+weight
    ch_d = ROOT.TChain(treename)
    ch_d.AddFile(filedict["Data" + untagged_suffix][channelcuts[channel]["Data"]]["path"])
    ch_d.Draw(todraw + ">> " + hist_data.GetName(), "("+ cut + " && "+ channelcuts[channel]["cut"] +")" + "*"+Mbtag_weight)

    ch_TT = ROOT.TChain(treename)
    ch_TT.AddFile(fileTT)
    ch_TT.Draw(todraw + ">> " + hist_TT.GetName(), finalcut)

    
    
    #hist_DY = ROOT.TH1F("untagged_DY_"+channel+"_%s"%(suffix), "untagged_DY_"+channel+"_%s"%(suffix), len(xbins)-1, xbins)
    hist_DY = hist_data.Clone()
    hist_DY.SetName("DY_"+channel+"_%s"%(suffix))
    hist_DY.Add(hist_TT, -1)
    print "Data-driven DY estimation, ch: ",channel," data ",hist_data.Integral(), " TT ",hist_TT.Integral()," DY ",hist_DY.Integral()
    #makeDYplots = False
    if makeDYplots:
        hs = ROOT.THStack("datadrivenDY_"+channel+"_"+suffix, "Data Driven Estimation for Drell-Yan")
        colors = [800-4, 820-3, 900-3, 860-3, 616-7, 432+2, 400+2]
        hist_TT.SetFillColor(colors[0])
        hist_DY.SetFillColor(colors[1])
        hist_data.SetMarkerColor(1)
        hist_data.SetMarkerStyle(20)
        hs.Add(hist_TT)
        hs.Add(hist_DY)
        legend = ROOT.TLegend(0.74,0.62,0.84,0.65+3*.05); 
	legend.SetTextSize(0.04); legend.SetTextFont(42)
        legend.SetBorderSize(0)
        legend.AddEntry(hist_data, "Data: %.1f"%hist_data.Integral(),"p")
        legend.AddEntry(hist_TT, "TT: %.1f"%hist_TT.Integral(),"f")
        legend.AddEntry(hist_DY, "DY = Data-TT","f")
 
        c1 =  ROOT.TCanvas()
        hs.Draw("hist")
        hist_data.Draw("epsame")
        legend.Draw("same")
        hs.GetHistogram().GetXaxis().SetTitle(xtitle)
        hs.GetHistogram().GetYaxis().SetTitle("Events")
	tex1 = ROOT.TLatex(0.17,0.8, channelcuts[channel]["latex"]+" channel, "+cut)
	tex1.SetNDC(); tex1.SetTextSize(.055)
        tex1.Draw("same")
        #plotdir = "DataDriven_DY_plots/"
        c1.SaveAs(plotname+"_"+channel+"_"+suffix+".pdf")
          
        
        
    hist_DY.SetDirectory(0)
    hist_data.SetDirectory(0)
    hist_TT.SetDirectory(0)

    return [hist_DY,hist_data, hist_TT]


def addTH1withError(hist1, hist2, c2=1.0):
    
    h_dummy = hist1.Clone()
    for bin in xrange(h_dummy.GetNbinsX()):
        value1 = hist1.GetBinContent(bin+1)
        err1 = hist1.GetBinError(bin+1)
        value2 = hist2.GetBinContent(bin+1)
        err2 = hist2.GetBinError(bin+1)
        h_dummy.SetBinContent(bin+1, value1+ value2*c2)
        h_dummy.SetBinError(bin+1, sqrt(err1*err1 + err2*err2*c2*c2))

    #h_dummy.Print("ALL")
    return h_dummy

def addTH2withError(hist1, hist2, c2=1.0):
    
    h_dummy = hist1.Clone()
    for bin in xrange(h_dummy.GetNbinsX()):
        for ybin in xrange(h_dummy.GetNbinsY()):
            value1 = hist1.GetBinContent(bin+1, ybin+1)
            err1 = hist1.GetBinError(bin+1, ybin+1)
            value2 = hist2.GetBinContent(bin+1, ybin+1)
            err2 = hist2.GetBinError(bin+1, ybin+1)
            h_dummy.SetBinContent(bin+1, ybin+1, value1+ value2*c2)
            h_dummy.SetBinError(bin+1, ybin+1, sqrt(err1*err1 + err2*err2*c2*c2))

    #h_dummy.Print("ALL")
    return h_dummy

def add_sytematic_statistic_error(hist):
    for bin in xrange(hist.GetNbinsX()):
        syserr = hist.GetBinError(bin+1)
        hist.SetBinError(bin+1, sqrt(syserr*syserr + abs(hist.GetBinContent(bin+1))))

    
###
def histForlimits1D(bgnames, mass, todraw, cut, xbins, xtitle, suffix, outfile, plotname):

    DYdatadriven = True
    if len(xbins) == 3:
        nbins = xbins[0]; xmin = xbins[1]; xmax =  xbins[2]
        xbins = []
        binwidth = (xmax-xmin)/nbins
        for i in range(0, nbins+1):
            xbins.append(xmin + i*binwidth)
        xbins = np.asarray(xbins)
    LouvainPlot = False

    signalmasses = [300, 400, 600, 900]
    signalcolor = [46,30,40,28, 20] 
    #print "nnout ",nnout
    treename = "Friends"
    #treename = "evtreeHME_nn"

    chlist = {}
    for shortname in full_local_samplelist.keys():
        for samplename in full_local_samplelist[shortname]:
            chlist[samplename] = ROOT.TChain(treename)
            chlist[samplename].AddFile(full_local_samplelist[shortname][samplename]['path'])
    
    signalname_short = 'RadionM%d'%mass; signalname_full =  "GluGluToRadionToHHTo2B2VTo2L2Nu_M-%d_narrow_13TeV-madgraph-v2"%mass
    filesignal = full_local_samplelist[signalname_short][signalname_full]["path"]
    ch_s =  ROOT.TChain(treename); ch_s.Add(filesignal)

    ch_s_toplot = []; signal_xsec = 10.0
    for smass in signalmasses:
        signalname_short = 'RadionM%d'%smass; signalname_full =  "GluGluToRadionToHHTo2B2VTo2L2Nu_M-%d_narrow_13TeV-madgraph-v2"%smass
        filesignal = full_local_samplelist[signalname_short][signalname_full]["path"]
        ch_s_toplot.append(ROOT.TChain(treename)); ch_s_toplot[-1].Add(filesignal)



    directory = plotname.replace(plotname.split("/")[-1],'')
    #colors = [628, 596, 820, 432, 418]
    colors = [800-4, 820-3, 900-3, 860-3, 616-7, 432+2, 400+2]
    nominalWeight = "sample_weight*event_reco_weight*{totallumi}".format(totallumi = TotalLumi)
    systematiclist = ["CMS_eff_b_heavy","CMS_eff_b_light","CMS_pu", "CMS_pdf", "CMS_eff_trigger","CMS_eff_e","CMS_eff_mu","CMS_iso_mu","QCDscale", "MC_statistical"]
    #systematiclist = ["QCDscale"]
    #systematiclist = ["CMS_eff_b_heavy"]
    #systematiclist = ["CMS_pdf"]
    
    sysplotter = bbWWPlotterSystematics(full_local_samplelist, nominalWeight, systematiclist, DYdatadriven)
    sysplotter.initialize1D(todraw, xtitle, xbins, cut)
    xsec_signal = 5.0#pb
    sysplotsCollector = {}
    for channel in channelcuts:
	allhists = []
        hists_datadriven = []
        allhists_v2 = {}
        sysplotsCollector[channel] = {}
        rfile = ROOT.TFile(outfile, "UPDATE")
        rfile.Close()
        legend = ROOT.TLegend(0.74,0.62,0.84,0.65+len(bgnames)*.04); 
	legend.SetTextSize(0.045); legend.SetTextFont(42)
        legend.SetBorderSize(0)
        legend2 = ROOT.TLegend(0.43,0.7,0.65,0.72+3*.05); 
	legend2.SetTextSize(0.035); legend2.SetTextFont(42)
        legend2.SetBorderSize(0)
        legend2.SetHeader("signal #sigma=%.1f pb"%(signal_xsec))
	BGSum = 0.0
        maxbgbin = 0.0
        #hist_data_fake = ROOT.TH1F("data_obs_"+channel+"_M%d_%s"%(mass, suffix), "data_obs_"+channel+"_M%d_%s"%(mass, suffix), len(xbins)-1, xbins)
	hist_data = ROOT.TH1F("data_obs_"+channel+"_%s"%(suffix), "data_obs_"+channel+"_%s"%( suffix), len(xbins)-1, xbins)
        ch_d = ROOT.TChain( treename )
        ch_d.AddFile(full_local_samplelist["Data"][channelcuts[channel]["Data"]]["path"])
        ch_d.Draw(todraw + ">> " + hist_data.GetName(), "("+ cut + " && "+ channelcuts[channel]["cut"] +")")
        hist_allsignal = []
        for im, smass in enumerate(signalmasses):
            hist_allsignal.append(ROOT.TH1F("signal_"+channel+"_M%d"%(smass), "signal_"+channel+"_M%d"%(smass), len(xbins)-1, xbins))
            weight = "(35.87*1000*%f/event_weight_sum)*final_total_weight"%(signal_xsec)
            ch_s_toplot[im].Draw(todraw + ">> " + hist_allsignal[im].GetName(), "("+ cut + " && "+ channelcuts[channel]["cut"] +")*"+weight)
            hist_allsignal[im].SetLineColor(signalcolor[im])
            hist_allsignal[im].SetLineWidth(2)
            legend2.AddEntry(hist_allsignal[im],"signal %d GeV"%smass, "l")




        sysplotter.runSystematics('RadionM%d'%mass, channel)
        hist_s = sysplotter.finalhist["nominal"]
        sysplotsCollector[channel]['RadionM%d'%mass] = sysplotter.channel_shortname_systematic_hist[channel]['RadionM%d'%mass].copy()
        
        #sysplotter.writeSystematicsToFile( directory)

	hist_s.SetLineColor(colors[-1]+10)
	hist_s.SetLineWidth(2)
        hist_data.SetMarkerStyle(20)
        hist_data.SetMarkerColor(1)
        hist_data.SetLineColor(1)
        #legend2.AddEntry(hist_data,"Data","p")
        #legend2.AddEntry(hist_s,"#splitline{Signal}{%d GeV, %d pb}"%(mass, xsec_signal),"l")


	hist_bg_all = ROOT.TH1F("bg_all_"+channel+"_%s"%(suffix), "bg_all_"+channel+"_%s"%( suffix), len(xbins)-1, xbins)
        for bin in xrange(hist_bg_all.GetNbinsX()):  
            hist_bg_all.SetBinContent(bin+1, 0.0)
            hist_bg_all.SetBinError(bin+1, 0.0)


	for i, key in enumerate(bgnames):
	    hist = ROOT.TH1F(key+"_"+channel+"_%s"%(suffix), key+"_"+channel+"_%s"%(suffix), len(xbins)-1, xbins)

            sysplotter.runSystematics(key, channel)
            if key == "DY" and DYdatadriven and (channel == "MuMu" or channel == "ElEl"):
                hists_datadriven.append(sysplotter.hist_data_untagged)
                hist = sysplotter.finalhist["nominal"]
                for key_datadriven in sysplotter.shortnames_datadriven:
                    sysplotsCollector[channel][key_datadriven] = sysplotter.channel_shortname_systematic_hist[channel][key_datadriven].copy()
                    #print "systematic hist for DY datadriven ",sysplotsCollector[channel][key_datadriven]
                    #hist = addTH1withError(hist, sysplotter., -1)
            else:
                sysplotsCollector[channel][key] = sysplotter.channel_shortname_systematic_hist[channel][key].copy()
                #hist = sysplotter.finalhist["nominal"]
                hist = sysplotsCollector[channel][key]["nominal"]
	    
            allhists.append(hist)
	    hist.SetFillColor(colors[i])
	    hist.SetLineColor(colors[i])
            
            #if key == "TT":
            #    maxbgbin = hist.GetBinContent(hist.GetMaximumBin())
	    BGSum = BGSum + hist.Integral()
            #hist_data.Add(hist)
            #hist_bg_all.Add(hist)
            hist_bg_all = addTH1withError(hist_bg_all, hist)
	    print "mass ",mass, " channel ", channel," bg ",key," rate ",hist.Integral()
            #hist_bg_all.Print("ALL")
            

        hist_bg_all.SetLineColor(0)
        hist_bg_all.SetFillStyle(3244)
        hist_bg_all.SetFillColor(14)
        #add_sytematic_statistic_error(hist_bg_all) ## MC statistic error is already considered
        #maxbgbin = hist_data.GetBinContent(hist_data.GetMaximumBin())
        #maxbgbin = hist_data.GetBinContent(hist_bg_all.GetMaximumBin())
        maxbgbin = hist_bg_all.GetBinContent(hist_bg_all.GetMaximumBin())
        maxsignalbin = hist_s.GetBinContent(hist_s.GetMaximumBin())
        #print "maxsignalbin ",maxsignalbin," maxbgbin ",maxbgbin, " max ",hist_bg_all.GetMaximumBin()


        rfile = ROOT.TFile(outfile, "UPDATE")
        hs = ROOT.THStack("allbg_"+channel+"_"+suffix, "  ")
	for i in range(len(allhists)):
            index = len(allhists)-i-1
	    hs.Add(allhists[index])
	    legend.AddEntry(allhists[i], bgnames[i], "f")
            #allhists[index].SetDirectory(rfile)
            #allhists[index].Write()
        for bgname in sysplotsCollector[channel].keys():
            sysplotsCollector[channel][bgname]["nominal"].SetDirectory(rfile)
            sysplotsCollector[channel][bgname]["nominal"].Write()
            for sys in sysplotsCollector[channel][bgname].keys():
                if sys == "nominal":
                    continue
                for plottype in ["up","down"]:
                    #print "channel ",channel, " bgname ", bgname," sys ",sys," ",plottype, " hist name ",sysplotsCollector[channel][bgname][sys][plottype].GetName()," ",sysplotsCollector[channel][bgname][sys][plottype].Print()

                    sysplotsCollector[channel][bgname][sys][plottype].SetDirectory(rfile)
                    sysplotsCollector[channel][bgname][sys][plottype].Write()

        if DYdatadriven:
            for h in hists_datadriven:
                h.SetDirectory(rfile)
                h.Write()


	hs.Write()
        if maxsignalbin>maxbgbin:
            hs.SetMaximum(maxsignalbin*1.4)
        else:
            hs.SetMaximum(maxbgbin*1.1)

        print "mass ",mass, " channel ", channel," rate: signal ",hist_s.Integral()," BG ",BGSum," data ",hist_data.Integral()," Data/MC ",hist_data.Integral()/BGSum," MC:S/sqrt(B) ",hist_s.Integral()/sqrt(BGSum)

        hist_s.SetDirectory(rfile)
        hist_bg_all.SetDirectory(rfile)
        hist_data.SetDirectory(rfile)
	hist_s.Write()
        hist_bg_all.Write()
	hist_data.Write()

        ROOT.gStyle.SetPadLeftMargin(0.16)
        c1 = ROOT.TCanvas("c", "canvas", 800, 800)
        c1.Clear()
        
        ### top plot
        #pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.0, 1, 1.0)
        #pad1.SetBottomMargin(.02)
        pad1.Draw()
        pad1.cd()
        hs.Draw("hist")
        hist_bg_all.Draw("e3same")
        #hist_s.Draw("samehist")
        for hist_sig in hist_allsignal:
            hist_sig.Draw("samehist")
        #hist_data.Draw("epsame")
	
        #tex1 = ROOT.TLatex(0.17,0.8, channelcuts[channel]["latex"]+" channel, "+nnout)
	tex1 = ROOT.TLatex(0.2,0.83, channelcuts[channel]["latex"]+" channel ")
	tex1.SetNDC(); tex1.SetTextSize(.045)
        #tex2 = ROOT.TLatex(0.19,0.6, "M_{jj}<75 GeV "+"  "*12+" 75 GeV <=M_{jj}<140 GeV"+ "  "*12+" M_{jj} >= 140 GeV")
	tex2 = ROOT.TLatex(0.19,0.55, "M_{jj}<75 GeV "+"  "*6+" 75 GeV <=M_{jj}<140 GeV"+ "  "*6+" M_{jj} >= 140 GeV")
	tex2.SetNDC(); tex2.SetTextSize(.035)
	tex1.Draw("same")
        #tex2.Draw("same")
        if ( LouvainPlot):
           tex2.Draw("same")
	#hs.GetHistogram().SetTitle(" #scale[1.4]{#font[61]{CMS}} Simulation Preliminary"+"  "*38+"35.87 fb^{-1} (13 TeV),2016")
        #hs.GetHistogram().SetTitle("")
	#hs.GetHistogram().SetTitleSize(.04)
	#hs.GetHistogram().SetTitleOffset(1.2)
        tex0 = ROOT.TLatex(0.1,0.92, " #scale[1.4]{#font[61]{CMS}} #font[42]{Work in progress}"+"  "*2+"#font[42]{35.87 fb^{-1} (13 TeV),2016}")
	tex0.SetNDC(); tex0.SetTextSize(.035)
	tex0.Draw("same")
        #hs.GetHistogram().GetXaxis().SetTitle("DNN output, M_{jj} bins")
        #hs.GetHistogram().GetXaxis().SetTitle(xtitle)
        ## remove x label for hs,  
        #xaxis = hs.GetHistogram().GetXaxis()
        #for i in xrange(hs.GetHistogram().GetNbinsX()):
        #    xaxis.SetBinLabel(i+1, "")
        hs.GetHistogram().GetYaxis().SetTitle("Events")
        hs.GetHistogram().GetYaxis().SetTitleSize(.05)
        hs.GetHistogram().GetYaxis().SetLabelFont(42)
        hs.GetHistogram().GetYaxis().SetLabelSize(.045)
        hs.GetHistogram().GetYaxis().SetTitleOffset(1.4)
        hs.GetHistogram().GetXaxis().SetTitleSize(.045)
        hs.GetHistogram().GetXaxis().SetLabelFont(42)
        hs.GetHistogram().GetXaxis().SetTitle(xtitle)
        legend.Draw("same")
        legend2.Draw("same")
	#tex.Draw("same")


        """
        c1.cd()
        c1.Update()


        ### bottom plot
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, .29)
        pad2.SetTopMargin(0.0)
        pad2.SetBottomMargin(.35)
        pad2.SetTicks(1,1)#draw x, y axis on both side (left right for y, and bottom up for x)
        pad2.SetGridy()
        pad2.Draw()
        pad2.cd()
        hratio_framework = hist_data.Clone(); hratio_framework.SetName("hratio_framework")
        hratio = hist_data.Clone(); hratio.SetName("ratio")
        hratio.SetMarkerStyle(20)
        hratio.SetMarkerColor(1)
        herrband = hist_data.Clone(); herrband.SetName("errband")
        herrband.SetFillStyle(3244)
        herrband.SetFillColor(14)
        herrband.SetMarkerColor(0)
        herrband.SetMarkerSize(0)
        hratio.Divide(hist_bg_all)
        for bin in xrange(hratio.GetNbinsX()):
            value_den = hist_bg_all.GetBinContent(bin+1)
            err_den = hist_bg_all.GetBinError(bin+1)
            ratio = hratio.GetBinContent(bin+1)
            herrband.SetBinContent(bin+1, 1.0)
            #print "bin ",bin+1, " BG hist entry ", value_den," error ", err_den," ratio ",ratio," error on ratio ",ratio*err_den/value_den
            if abs(value_den) > 0.0:
                herrband.SetBinError(bin+1, ratio*err_den/value_den)
            else:
                herrband.SetBinError(bin+1, 0.0)
            hratio_framework.SetBinContent(bin+1, -1.)
        deltaY = 0.39## rather than use 0.4 to avoid
        hratio_framework.Draw()
        herrband.Draw("e3same")
        hratio.Draw("psame")

        herrband.SetStats(0)
        hratio.SetStats(0)

        hratio_framework.SetTitle("")
        hratio_framework.SetMaximum(1.0 + deltaY)
        hratio_framework.SetMinimum(1.0 - deltaY)
        hratio_framework.GetXaxis().SetTitle(xtitle)
        hratio_framework.GetXaxis().SetTitleSize(35)
        hratio_framework.GetXaxis().SetTitleFont(43)
        hratio_framework.GetXaxis().SetTitleOffset(3.0)
        hratio_framework.GetXaxis().SetLabelSize(30)
        hratio_framework.GetXaxis().SetLabelFont(43)#Absolute font size in pixel (precision 3)
        hratio_framework.GetYaxis().SetTitle("Data/MC")
        hratio_framework.GetYaxis().SetNdivisions(505)
        hratio_framework.GetYaxis().CenterTitle()
        hratio_framework.GetYaxis().SetTitleSize(25)
        hratio_framework.GetYaxis().SetTitleFont(43)
        hratio_framework.GetYaxis().SetTitleOffset(1.3)
        hratio_framework.GetYaxis().SetLabelSize(25)
        hratio_framework.GetYaxis().SetLabelFont(43)#Absolute font size in pixel (precision 3)

        #tex_pad2 = ROOT.TLatex(0.2,0.35, "Maximum #frac{S}{#sqrt{B}} = %.1f @ %.3f"%(bestS, bestWP))
        #tex_pad2.SetNDC()
        #tex_pad2.SetTextSize(.035)
        #tex_pad2.Draw("same")

        """

        c1.SaveAs(plotname+"_Radion_"+channel+"_"+suffix+".C")
        #c1.SaveAs(plotname+"_Radion_"+channel+".png")
        c1.SaveAs(plotname+"_Radion_"+channel+"_"+suffix+".pdf")
        rfile.Close()

        print "done with histForlimits1D"
	

def plotDNNoutput(masspoints, nnout, outdir):
    if not os.path.isdir(outdir):
        os.system("mkdir -p "+outdir)
    bgnames = ["TT","DY","sT","VV", "ttV"]
    nnbins = []
    for x in range(3, 75-6+1):
        nnbins.append(x*1.0/25.0)
    #print "nnbins ",nnbins
    nnbins_x = np.asarray(nnbins)
    xtitle = "DNN output, Mjj bin"
    for mass in masspoints:
        outfile = os.path.join(outdir, "Hhh_FinalBGYield_xsec1pb_NN_%s_SignalM%d.root"%( nnout, mass))
        plotname = os.path.join(outdir, "Hhh_FinalBGYield_xsec1pb_NN_%s"%nnout)
        tfile = ROOT.TFile(outfile, "RECREATE")
        tfile.Close()

        #cut = "({nnout}_M{mass}>3.0/25 && hme_h2mass_reco>=250)".format(nnout = nnout, mass = mass)
        cut = "({nnout}_M{mass}>3.0/25)".format(nnout = nnout, mass = mass) ##no HME cut for Louvain case
        todraw = "(({nnout}_M{mass})*(jj_M<75 && {nnout}_M{mass}>3.0/25)+(jj_M>=75 && jj_M<140 && {nnout}_M{mass}>3.0/25)*({nnout}_M{mass}+1-3.0/25)+(jj_M>=140 && {nnout}_M{mass}>3.0/25)*({nnout}_M{mass}+2-6.0/25))".format(nnout = nnout, mass=mass)
        suffix = "M%d"%(mass)
        histForlimits1D(bgnames, mass, todraw, cut, nnbins_x, xtitle, suffix, outfile, plotname)


def makeBackgroundshist(masspoints, variable, nbins, xtitle, outdir):

    def makeDYEstimationplots():
        for channel in ["MuMu","ElEl"]:
	    plotdir = "dataDriven_DYestimation/"
            plotname1 = os.path.join(plotdir, "Kinematics_%s"%variable+"_llMLT76")
            makeDYplots = True
            #plotname2 = os.path.join(plotdir, "Kinematics_%s"%variable+"_llMGT76")
            DrellYanDataDriven(channel, full_local_samplelist, variable, "ll_M<76", nbins, xtitle, "v0", makeDYplots, plotname1)
            #DrellYanDataDriven(channel, untagged_samplelist, variable, "ll_M>76", nbins, xtitle, "v0", plotname2)
    #bgnames = ["TT","DY","sT","Wjet","VV","ttV"]
    #makeDYEstimationplots()
    #pass
    #bgnames = ["TT","DY","sT","VV", "Wjet","ttV"]
    bgnames = ["TT","DY","sT","VV", "ttV"]
    #bgnames = ["TT","DY", "Wjet","ttV"]
    #bgnames = ["TT"]
    plotname = os.path.join(outdir, "Kinematics_%s"%variable)
    ###create tfile
    todraw = variable
    for mass in masspoints:
        outfile = os.path.join(outdir, "Backgrounds_signalM%d_allinputs.root"%mass)
        tfile = ROOT.TFile(outfile, "RECREATE")
        tfile.Close()
        cut = "ll_M<76"
        #todraw = "(({nnout}_M{mass}-3.0/25)*(jj_M<75 && {nnout}_M{mass}>3.0/25)+(jj_M>=75 && jj_M<140 && {nnout}_M{mass}>3.0/25)*({nnout}_M{mass}+1-6.0/25)+(jj_M>=140 && {nnout}_M{mass}>3.0/25)*({nnout}_M{mass}+2-9.0/25))".format(nnout = nnout, mass=mass)
        suffix = 'M%d'%mass
        histForlimits1D(bgnames, mass, todraw, cut, nbins, xtitle, suffix, outfile, plotname)




    
###
def histForlimits2D(bgnames, mass, todraw, cut, xbins, xtitle, ybins, ytitle, suffix, outfile, plotname):

    DYdatadriven = True
    if len(xbins) == 3:
        nbins = xbins[0]; xmin = xbins[1]; xmax =  xbins[2]
        xbins = []
        binwidth = (xmax-xmin)/nbins
        for i in range(0, nbins+1):
            xbins.append(xmin + i*binwidth)
        xbins = np.asarray(xbins)

    if len(ybins) == 3:
        nbins = ybins[0]; ymin = ybins[1]; ymax =  ybins[2]
        ybins = []
        binwidth = (ymax-ymin)*1.0/nbins
        for i in range(0, nbins+1):
            ybins.append(ymin + i*binwidth)
        ybins = np.asarray(ybins)


    LouvainPlot = False

    #print "nnout ",nnout
    treename = "Friends"
    #treename = "evtreeHME_nn"

    chlist = {}
    for shortname in full_local_samplelist.keys():
        for samplename in full_local_samplelist[shortname]:
            chlist[samplename] = ROOT.TChain(treename)
            chlist[samplename].AddFile(full_local_samplelist[shortname][samplename]['path'])
    
    signalname_short = 'RadionM%d'%mass; signalname_full =  "GluGluToRadionToHHTo2B2VTo2L2Nu_M-%d_narrow_13TeV-madgraph-v2"%mass
    filesignal = full_local_samplelist[signalname_short][signalname_full]["path"]
    ch_s =  ROOT.TChain(treename); ch_s.Add(filesignal)


    directory = plotname.replace(plotname.split("/")[-1],'')
    #colors = [628, 596, 820, 432, 418]
    #colors = [800-4, 820-3, 900-3, 860-3, 616-7, 432+2, 400+2]
    nominalWeight = "sample_weight*event_reco_weight*{totallumi}".format(totallumi = TotalLumi)
    systematiclist = ["CMS_eff_b_heavy","CMS_eff_b_light","CMS_pu", "CMS_pdf", "CMS_eff_trigger","CMS_eff_e","CMS_eff_mu","CMS_iso_mu","QCDscale", "MC_statistical"]
    #systematiclist = ["QCDscale"]
    #systematiclist = ["CMS_eff_b_heavy"]
    #systematiclist = ["CMS_pdf"]
    
    sysplotter = bbWWPlotterSystematics(full_local_samplelist, nominalWeight, systematiclist, DYdatadriven)
    sysplotter.initialize2D(todraw, xtitle, xbins, ytitle, ybins, cut)
    xsec_signal = 1.0e-3#pb
    sysplotsCollector = {}
    for channel in channelcuts:
	allhists = []
        hists_datadriven = []
        allhists_v2 = {}
        sysplotsCollector[channel] = {}
        rfile = ROOT.TFile(outfile, "UPDATE")
        rfile.Close()
	BGSum = 0.0
        maxbgbin = 0.0
	hist_data = ROOT.TH2F("data_obs_"+channel+"_%s"%(suffix), "data_obs_"+channel+"_%s"%( suffix), len(xbins)-1, xbins, len(ybins)-1, ybins)
        ch_d = ROOT.TChain( treename )
        ch_d.AddFile(full_local_samplelist["Data"][channelcuts[channel]["Data"]]["path"])
        ch_d.Draw(todraw + ">> " + hist_data.GetName(), "("+ cut + " && "+ channelcuts[channel]["cut"] +")")


        sysplotter.runSystematics2D('RadionM%d'%mass, channel)
        hist_s = sysplotter.finalhist["nominal"]
        sysplotsCollector[channel]['RadionM%d'%mass] = sysplotter.channel_shortname_systematic_hist[channel]['RadionM%d'%mass].copy()
        
        #sysplotter.writeSystematicsToFile( directory)



	hist_bg_all = ROOT.TH2F("bg_all_"+channel+"_%s"%(suffix), "bg_all_"+channel+"_%s"%( suffix), len(xbins)-1, xbins, len(ybins)-1, ybins)
        for bin in xrange(hist_bg_all.GetNbinsX()):  
            for ybin in xrange(hist_bg_all.GetNbinsY()):  
                hist_bg_all.SetBinContent(bin+1, ybin+1, 0.0)
                hist_bg_all.SetBinError(bin+1, ybin+1, 0.0)


	for i, key in enumerate(bgnames):
	    hist = ROOT.TH2F(key+"_"+channel+"_%s"%(suffix), key+"_"+channel+"_%s"%(suffix), len(xbins)-1, xbins, len(ybins)-1, ybins)

            sysplotter.runSystematics2D(key, channel)
            if key == "DY" and DYdatadriven and (channel == "MuMu" or channel == "ElEl"):
                hists_datadriven.append(sysplotter.hist_data_untagged)
                hist = sysplotter.finalhist["nominal"]
                for key_datadriven in sysplotter.shortnames_datadriven:
                    sysplotsCollector[channel][key_datadriven] = sysplotter.channel_shortname_systematic_hist[channel][key_datadriven].copy()
                    #print "systematic hist for DY datadriven ",sysplotsCollector[channel][key_datadriven]
                    #hist = addTH2withError(hist, sysplotter., -1)
            else:
                sysplotsCollector[channel][key] = sysplotter.channel_shortname_systematic_hist[channel][key].copy()
                #hist = sysplotter.finalhist["nominal"]
                hist = sysplotsCollector[channel][key]["nominal"]
	    
            allhists.append(hist)
            
            #if key == "TT":
            #    maxbgbin = hist.GetBinContent(hist.GetMaximumBin())
	    BGSum = BGSum + hist.Integral()
            #hist_data.Add(hist)
            #hist_bg_all.Add(hist)
            hist_bg_all = addTH2withError(hist_bg_all, hist)
	    print "mass ",mass, " channel ", channel," bg ",key," rate ",hist.Integral()
            #hist_bg_all.Print("ALL")
            

        #add_sytematic_statistic_error(hist_bg_all) ## MC statistic error is already considered
        maxbgbin = hist_data.GetBinContent(hist_data.GetMaximumBin())
        maxsignalbin = hist_s.GetBinContent(hist_s.GetMaximumBin())


        rfile = ROOT.TFile(outfile, "UPDATE")
        for bgname in sysplotsCollector[channel].keys():
            sysplotsCollector[channel][bgname]["nominal"].SetDirectory(rfile)
            sysplotsCollector[channel][bgname]["nominal"].Write()
            for sys in sysplotsCollector[channel][bgname].keys():
                if sys == "nominal":
                    continue
                for plottype in ["up","down"]:
                    #print "channel ",channel, " bgname ", bgname," sys ",sys," ",plottype, " hist name ",sysplotsCollector[channel][bgname][sys][plottype].GetName()," ",sysplotsCollector[channel][bgname][sys][plottype].Print()

                    sysplotsCollector[channel][bgname][sys][plottype].SetDirectory(rfile)
                    sysplotsCollector[channel][bgname][sys][plottype].Write()

        if DYdatadriven:
            for h in hists_datadriven:
                h.SetDirectory(rfile)
                h.Write()



        print "mass ",mass, " channel ", channel," rate: signal ",hist_s.Integral()," BG ",BGSum," data ",hist_data.Integral()," Data/MC ",hist_data.Integral()/BGSum," MC:S/sqrt(B) ",hist_s.Integral()/sqrt(BGSum)

        hist_s.SetDirectory(rfile)
        hist_bg_all.SetDirectory(rfile)
        hist_data.SetDirectory(rfile)
	hist_s.Write()
        hist_bg_all.Write()
	hist_data.Write()

        ROOT.gStyle.SetPadLeftMargin(0.18)
        ROOT.gStyle.SetPadRightMargin(0.15)
        c1 = ROOT.TCanvas("c", "canvas", 800, 800)
        c1.Clear()
        c1.Divide(2, 2)
        

        c1.cd(1)
        hist_s.SetTitle("Signal, M=%d GeV, 2016"%mass)
	hist_s.Draw("colz")
        hist_s.GetXaxis().SetTitle(xtitle)
        hist_s.GetYaxis().SetTitle(ytitle)
        hist_s.GetXaxis().CenterTitle()
        hist_s.GetYaxis().CenterTitle()
        hist_s.GetXaxis().SetTitleSize(.05)
        hist_s.GetXaxis().SetLabelFont(42)
        hist_s.GetYaxis().SetTitleSize(.05)
        hist_s.GetYaxis().SetLabelFont(42)
        hist_s.GetYaxis().SetLabelSize(.045)
        hist_s.GetYaxis().SetTitleOffset(1.3)
        c1.Update()
        c1.cd(2)
        hist_bg_all.SetTitle("SM backgrounds, 2016")
        hist_bg_all.Draw("colz")
        hist_bg_all.GetXaxis().SetTitle(xtitle)
        hist_bg_all.GetYaxis().SetTitle(ytitle)
        hist_bg_all.GetXaxis().CenterTitle()
        hist_bg_all.GetYaxis().CenterTitle()
        hist_bg_all.GetXaxis().SetTitleSize(.05)
        hist_bg_all.GetXaxis().SetLabelFont(42)
        hist_bg_all.GetXaxis().SetLabelSize(.045)
        hist_bg_all.GetYaxis().SetTitleSize(.05)
        hist_bg_all.GetYaxis().SetLabelFont(42)
        hist_bg_all.GetYaxis().SetLabelSize(.045)
        hist_bg_all.GetYaxis().SetTitleOffset(1.3)
        c1.Update()
        c1.cd(3)
        hist_data.SetTitle("Data, 2016")
        hist_data.Draw("colz")
        hist_data.GetXaxis().SetTitle(xtitle)
        hist_data.GetYaxis().SetTitle(ytitle)
        hist_data.GetXaxis().CenterTitle()
        hist_data.GetYaxis().CenterTitle()
        hist_data.GetXaxis().SetTitleSize(.05)
        hist_data.GetXaxis().SetLabelFont(42)
        hist_data.GetXaxis().SetLabelSize(.045)
        hist_data.GetYaxis().SetTitleSize(.05)
        hist_data.GetYaxis().SetLabelFont(42)
        hist_data.GetYaxis().SetLabelSize(.045)
        hist_data.GetYaxis().SetTitleOffset(1.3)
        c1.Update()
        c1.cd(4)
        histRatio = hist_data.Clone()
        histRatio.Divide(hist_bg_all)
        histRatio.SetTitle("Data/MC, 2016")
        histRatio.GetXaxis().CenterTitle()
        histRatio.GetYaxis().CenterTitle()
        histRatio.GetXaxis().SetTitleSize(.05)
        histRatio.GetXaxis().SetLabelFont(42)
        histRatio.GetXaxis().SetLabelSize(.045)
        histRatio.GetYaxis().SetTitleSize(.05)
        histRatio.GetYaxis().SetLabelFont(42)
        histRatio.GetYaxis().SetLabelSize(.045)
        histRatio.GetYaxis().SetTitleOffset(1.3)
        histRatio.Draw("colz")
        maxerr  = 0.2
        histRatio.SetMaximum(1.+ maxerr)
        histRatio.SetMinimum(1.- maxerr)
        c1.Update()
        c1.cd()

        c1.SaveAs(plotname+"_Radion_"+channel+"_"+suffix+".C")
        #c1.SaveAs(plotname+"_Radion_"+channel+".png")
        c1.SaveAs(plotname+"_Radion_"+channel+"_"+suffix+".pdf")
        rfile.Close()

        print "done with histForlimits1D"



def generateHMEbins(mass):
    lowM = 250.0; highM = 1200.0
    if highM > 2*mass:
       highM = 2*mass+100
    if lowM < mass/2.0:
       lowM = mass/2.0
    if mass == 750:
        lowM = 400
    xbins = [lowM]
    x = lowM
    step = 25.0 
    if mass >= 400 and mass < 500:
        step = 30.0
    elif mass>=500 and mass < 700:
        step = 40.0
    elif mass>=700:
        step = 60.0
    gap1 = 40 + mass*.1 #120.0
    gap2 = 50 + mass*.2 #200.0
    gap3 = 60 + mass*.3 #
    while x <= highM-150.0:
        if abs(x-mass) > mass:
            x = x+100.0
            xbins.append(x)
        elif abs(x-mass)<= gap2 and abs(x-mass)>gap1:
            x = x+step*1.6
            xbins.append(x)
        elif abs(x-mass)> gap3:
            x = x+step*3
            xbins.append(x)
        #elif abs(x-mass)> 200.:
        #    x = x+step*2.5
        #    xbins.append(x)
        else:
            x = x+step
            xbins.append(x)
    if mass<800:
        xbins = xbins[:-1]
    xbins.append(highM)
    print "Benchmark ", mass ," HME mass bins ",xbins

    if mass == 400:
        xbins = [250.0, 310.0, 350.0, 385.0, 415.0, 445.0, 475.0, 520.0, 570.0, 650.0, 900.0]
    return np.asarray(xbins)
            
	

def plotDNNoutputvsHME(masspoints, nnout, outdir):
    if not os.path.isdir(outdir):
        os.system("mkdir -p "+outdir)
    bgnames = ["TT","DY","sT","VV", "ttV"]


    xtitle = "DNN output, Mjj bin"
    ytitle = "HME"
    for mass in masspoints:
        hmebins = generateHMEbins(mass)
        nnstep = 0.05
        nncut = 0.35
        if mass>750:
            nncut = 0.35 + (mass-700)*.001 
        elif mass == 750:
            nncut = 0.45
        nn_min = nncut
        nn_max = 3.0 - nncut*2
        nnbins_x = np.arange(nn_min, nn_max+nnstep, nnstep)
        print "nnbins_x ",nnbins_x
        cut = "({nnout}_M{mass} > {nncut}) && hme_h2mass_reco>=250".format(nnout = nnout, mass = mass, nncut = nncut)
        xtodraw  = "(({nnout}_M{mass})*(jj_M<75) + (jj_M>=75 && jj_M<140)*({nnout}_M{mass}+1-1*{nncut})+(jj_M>=140)*({nnout}_M{mass}+2- 2*{nncut}))".format(nnout = nnout, mass=mass, nncut = nncut)
        outfile = os.path.join(outdir, "Hhh_FinalBGYield_xsec1pb_NNvsHME_%s_SignalM%d.root"%( nnout, mass))
        plotname = os.path.join(outdir, "Hhh_FinalBGYield_xsec1pb_NNvsHME_%s"%nnout)
        tfile = ROOT.TFile(outfile, "RECREATE")
        tfile.Close()

        todraw = "hme_h2mass_reco : "+xtodraw
        suffix = "M%d"%(mass)
        histForlimits2D(bgnames, mass, todraw, cut, nnbins_x, xtitle, hmebins, ytitle, suffix, outfile, plotname)




variablesdir = "HHNtuple_20180703_variablehist_addSys_HHworkshop/"
os.system("mkdir -p "+variablesdir)
varibales = ['jj_pt', 'll_pt', 'll_M', 'll_DR_l_l', 'jj_DR_j_j', 'llmetjj_DPhi_ll_jj', 'llmetjj_minDR_l_j', 'llmetjj_MTformula','mt2', 'jj_M','hme_h2mass_reco']
mcnames = ["TT","DY","sT","Wjet","VV","ttV"]
masspoints = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900]
#variables = ['lep1_pt']
#makeBackgroundshist(output_folder, [400], 'llmetjj_MTformula', [50, 0.0, 500.0],"MT", variablesdir)
def plotallkinematics():
    #output_folder = "/Users/taohuang/Documents/DiHiggs/20180316_NanoAOD/HHNtuple_20180328_fixedleptonDZeff"
    #print "Ntuple folder ",output_folder
    
    #makeBackgroundshist(masspoints, 'll_M', [50, 12.0, 76.0], "M_{ll}", variablesdir)
    #makeBackgroundshist([400], 'lep1_pt', [60, 10.0, 200], "lep1 p_{T}", variablesdir)
    #makeBackgroundshist([400], 'lep2_pt', [60, 10.0, 200], "lep2 p_{T}", variablesdir)
    #makeBackgroundshist([400], 'lep1_eta', [60, -2.4, 2.4], "lep1 #eta", variablesdir)
    #makeBackgroundshist([400], 'lep2_eta', [60, -2.4, 2.4], "lep2 #eta", variablesdir)
    #makeBackgroundshist([400], 'jet1_pt', [70, 20.0, 300], "jet1 p_{T}", variablesdir)
    #makeBackgroundshist([400], 'jet2_pt', [70, 20.0, 300], "jet2 p_{T}", variablesdir)
    #makeBackgroundshist([400], 'jet1_eta', [70, -2.5, 2.5], "jet1 #eta", variablesdir)
    #makeBackgroundshist([400], 'jet2_eta', [70, -2.5, 2.5], "jet2 #eta", variablesdir)
    #makeBackgroundshist([400], 'met_pt', [50, 0.0, 500.0],"MET p_{T}", variablesdir)
    ##makeBackgroundshist([400], 'met_phi', [60, -3.2, 3.20],"MET #phi", variablesdir)
    #makeBackgroundshist([400], 'll_M', [50, 12.0, 76.0], "M_{ll}", variablesdir)
    #makeBackgroundshist([400], 'll_DR_l_l', [50, .0, 6.0], "#DeltaR_{ll}", variablesdir)
    #makeBackgroundshist([400], 'jj_M', [50, 0.0, 400.0], "M_{jj}",variablesdir)
    #makeBackgroundshist([400], 'jj_DR_j_j', [50, .0, 6.0], "#DeltaR_{jj}",variablesdir)
    #makeBackgroundshist([400], 'abs(llmetjj_DPhi_ll_jj)', [24, .0, 3.1415926],"#Delta#phi(ll,jj)", variablesdir)
    #makeBackgroundshist([400], 'll_pt', [50, 0.0, 450.0], "Dilepton p_{T}", variablesdir)
    #makeBackgroundshist([400], 'jj_pt', [50, 0.0, 450.0], "Dijet p_{T}", variablesdir)
    #makeBackgroundshist([400], 'llmetjj_minDR_l_j', [50, .0, 5.0], "#DeltaR_{l,j}", variablesdir)
    #makeBackgroundshist([400], 'llmetjj_MTformula', [50, 0.0, 500.0],"MT", variablesdir)
    #makeBackgroundshist([400], 'mt2', [50, 0.0, 300.0],"MT2", variablesdir)
    #makeBackgroundshist([400], 'hme_h2mass_reco', [50, 250.0, 700.0],"HME", variablesdir)
    makeBackgroundshist([400], 'hme_h2mass_reco', [80, 250.0, 1200.0],"HME", variablesdir)
    #makeBackgroundshist([400], 'nnout_MTonly_M400', [20, 0.20, 1.0],"DNN output", variablesdir)


for mass in [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900]:
    generateHMEbins(mass)

masspoints = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900]
#masspoints = [300, 400, 600, 900]
masspoints = [300]
#masspoints = [900]
plotallkinematics()
nnout_list = ["nnout_MTonly","nnout_MTandMT2","nnout_MTandMT2_MJJ","nnout_MTandMT2_HME"]
#outdir = "HHbbWW_20180625_NNoutput_MjjCR_test/"
outdir = "HHbbWW_20180820_NNoutput_MjjCR_HHworkshop/"
#plotDNNoutput(masspoints, "nnout_MTonly",outdir)
#outdir = "HHbbWW_20180802_NNoutputVsHME_MjjCR_M400/"
#plotDNNoutputvsHME(masspoints, "nnout_MTonly", outdir)
bgnames = ["TT","DY","sT","Wjet","VV","ttV"]
#bgnames = ["TT","DY","sT","VV","ttV"]
#outcutflowdir = "HHNtuple_20180412_cutflows_newTT/"
outcutflowdir = "HHNtuple_20180518_addSys_cutflows/"
#os.system("mkdir -p "+outcutflowdir)
#plotCutflowHist(outcutflowdir, "TT")
#plotCutflowHist_data(outcutflowdir)
for mass in masspoints:
    mcnames.append("RadionM%d"%mass)
#plotCutflowHist_allMC(outcutflowdir, bgnames)
#runallCutflowhist(outcutflowdir, mcnames)
