import os
import numpy as np
import ROOT
from math import sqrt

import sys
sys.argv.append( '-b' )
sys.argv.append( '-q' )
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)


TotalLumi = 36.8 #fb-1
#TotalLumi = 35.920
#extraweight is used for compensating the "HLT" safe ID cut, for Louvain ntuples
channelcuts = {"MuMu":{"cut":"isMuMu","extraweight": 1.0,"Data":"DoubleMuon", "latex":"#mu#mu"},
        "ElMu":{"cut":"(isElMu)", "extraweight": 1.0,"Data":"MuonEG", "latex":"e#mu"},
        "MuEl":{"cut":"(isMuEl)", "extraweight": 1.0,"Data":"MuonEG", "latex":"#mue"},
        "ElEl":{"cut":"isElEl","extraweight":1.0,"Data":"DoubleEG", "latex":"ee"},
		}
#cutflows = ["All","Trigger","online-offline matching","dilepton PtEta","dilepton IP","dilepton ID","dilepton Iso","HLT Safe ID","nlepton>=3 veto","M_{ll}>12","NJets>=2","dijet PtEta","DR_j_l > 0.3","dijet btagging","M_{ll}<76"]
#colors = [623, 593, 820, 797, 432, 418]
colors = [800-4, 593, 820, 616-7, 432, 418]
bghistnamedict = {"TT":"ttbar","ttV":"ttV", "singleTop":"SingleTop","SMH":"SMHiggs","VV":"VV","Drell-Yan":"dy_mc"}
def plotshapes(file1, treename1, todraw1, cut1, file2, treename2, todraw2, cut2, xbins, xtitle, legs, plotname):
    
    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2]
    ch = plotname.split("_")[-1]
    

    ch1 = ROOT.TChain(treename1);  ch1.Add(file1)
    ch2 = ROOT.TChain(treename2);  ch2.Add(file2)
    print "file1 ",file1," todraw1 ",todraw1, " cut1 ",cut1
    print "file1 ",file2," todraw1 ",todraw2, " cut1 ",cut2
    h1 = ROOT.TH1F("hist1","hist1",xbins[0], xbins[1], xbins[2])
    h2 = ROOT.TH1F("hist2","hist2",xbins[0], xbins[1], xbins[2])
    ##weight include in cuts
    ch1.Draw(todraw1+">> hist1",cut1)
    ch2.Draw(todraw2+">> hist2",cut2)

    hs = ROOT.THStack("Shapecomparison", "  ")
    legend = ROOT.TLegend(0.65,0.7,0.88,0.7+2*.05); 
    legend.SetTextSize(0.04); legend.SetTextFont(42); legend.SetBorderSize(0)
    h1.SetLineColor(colors[0])
    h1.SetLineWidth(2)
    h2.SetLineColor(colors[1])
    h2.SetLineWidth(2)

    hs.Add(h1)
    hs.Add(h2)

    legend.AddEntry(h1, legs[0]+", Mean: %.2f"%h1.GetMean(), "l")
    legend.AddEntry(h2, legs[1]+", Mean: %.2f"%h2.GetMean(), "l")
    c1 = ROOT.TCanvas("c1","c1",800,600)
    c1.Clear()
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(.0)
    pad1.Draw()
    pad1.cd()
    hs.Draw("nostackhist")
    
    legend.Draw("same")
    tex1 = ROOT.TLatex(0.17,0.8, channelcuts[ch]["latex"]+ " channel")
    tex1.SetNDC(); tex1.SetTextSize(.045)
    tex1.Draw("same")
    tex2 = ROOT.TLatex(0.5,0.85, "Louvain: %d ,Mine: %d "%(h1.GetEntries(), h2.GetEntries()))
    tex2.SetNDC(); tex2.SetTextSize(.04)
    tex2.SetTextFont(42)
    tex2.Draw("same")
    c1.cd()
    c1.Update()


    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, .29)
    pad2.SetTopMargin(0.)
    pad2.SetBottomMargin(.35)
    pad2.SetTicks(1,1)#draw x, y axis on both side (left right for y, and bottom up for x)
    pad2.SetGridy()
    pad2.Draw()
    pad2.cd()

    hist_ratio = h1.Clone()
    hist_ratio.Divide(h2)

    hist_ratio.Draw("ep")
    hist_ratio.SetStats(0)
    hist_ratio.SetTitle("")
    maxratiobin = hist_ratio.GetBinContent(hist_ratio.GetMaximumBin())
    minratiobin = hist_ratio.GetBinContent(hist_ratio.GetMinimumBin())
    deltaY = abs(maxratiobin -1.0)
    if deltaY < abs(minratiobin -1.0):
        deltaY = abs(minratiobin -1.0)
    hist_ratio.GetYaxis().SetRangeUser(1-deltaY*1.3, 1+deltaY*1.3) 

    hist_ratio.GetXaxis().SetTitle(xtitle)
    hist_ratio.GetXaxis().SetTitleSize(20)
    hist_ratio.GetXaxis().SetTitleFont(43)
    hist_ratio.GetXaxis().SetTitleOffset(3.0)
    hist_ratio.GetXaxis().SetLabelSize(15)
    hist_ratio.GetXaxis().SetLabelFont(43)#Absolute font size in pixel (precision 3)
    hist_ratio.GetYaxis().SetTitle("Louvain/Mine")
    hist_ratio.GetYaxis().SetNdivisions(505)
    hist_ratio.GetYaxis().CenterTitle()
    hist_ratio.GetYaxis().SetTitleSize(20)
    hist_ratio.GetYaxis().SetTitleFont(43)
    hist_ratio.GetYaxis().SetTitleOffset(.9)
    hist_ratio.GetYaxis().SetLabelSize(15)
    hist_ratio.GetYaxis().SetLabelFont(43)#Absolute font size in pixel (precision 3)


    c1.SaveAs(plotname+".pdf")
    c1.SaveAs(plotname+".C")
    


def get_event_weight_sum_xsec(file1, isLouvain):
    tfile = ROOT.TFile(file1, "READ")
    event_weight_sum = 1.0; xsec = 1.0
    if isLouvain:
        a = tfile.Get("event_weight_sum")
        event_weight_sum = a.GetVal()
        b = tfile.Get("cross_section")
        xsec = b.GetVal()
    else:
        h_cutflow = tfile.Get("h_cutflow")
        event_weight_sum = h_cutflow.GetBinContent(1)
        b = tfile.Get("cross_section")
        xsec = b.GetVal()
        
    return (event_weight_sum, xsec)
        
        
def compareLouvainandMine(file1, file2, outputdir):

    treename1 = "t"
    treename2 = "Friends"
    todraw1_list = ["trigeff","llidiso","jjbtag_heavy*jjbtag_light", "event_pu_weight"]
    todraw2_list = ["lltrigger_weight","lliso_weight*llid_weight","event_btag_weight", "event_pu_weight"]
    llcut = "(ll_M<76)"

    event_weight_sum1, xsec1 = get_event_weight_sum_xsec(file1, True)
    event_weight_sum2, xsec2 = get_event_weight_sum_xsec(file2, False)
    #print "Louvain ",event_weight_sum1, xsec1," Mine ",event_weight_sum2, xsec2

    totalweight = "event_weight * trigeff * jjbtag_heavy * jjbtag_light * llidiso"# * pu"
    #weight1 = totalweight*"*1000.0*{totLumi}*{xsec}/{weight_sum}".format(totLumi = TotalLumi, xsec = xsec1, weight_sum = event_weight_sum1)
    #weight2 = totalweight*"sample_weight*1000.0*{totLumi}*{xsec}/{weight_sum}".format(totLumi = TotalLumi, xsec = xsec2, weight_sum = event_weight_sum2)
    weight1 = "1";weight2="1"

    
    legs = ["Louvain","Mine"]
    for ch in channelcuts:
        chcut = channelcuts[ch]["cut"]
        cut1 = "("+chcut + " && " + llcut+")*"+ weight1
        cut2 = "("+chcut + " && " + llcut+")*"+ weight2
        #print "cut1 ",cut1, " cut2 ",cut2
        for i in range(len(todraw1_list)):
            xbins = [60, .7, 1.30]
            todraw1 = todraw1_list[i]
            todraw2 = todraw2_list[i]
            xtitle = todraw2
            name = todraw2.replace("*","_")
            if todraw2 == "event_pu_weight":
                xbins = [60, 0.0, 3.0]
            plotname = os.path.join(outputdir, name+"_"+ch)
            plotshapes(file1, treename1, todraw1, cut1, file2, treename2, todraw2, cut2, xbins, xtitle, legs, plotname)


def compareKinematic(file1,file2, todraw, xbins, xtitle, variablesdir):
    treename1 = "t"
    treename2 = "Friends"
    llcut = "(ll_M<76)"

    event_weight_sum1, xsec1 = get_event_weight_sum_xsec(file1, True)
    event_weight_sum2, xsec2 = get_event_weight_sum_xsec(file2, False)
    print "Louvain ",event_weight_sum1, xsec1," Mine ",event_weight_sum2, xsec2

    totalweight = "event_weight * trigeff * jjbtag_heavy * jjbtag_light * llidiso"# * pu"
    #weight1 = totalweight*"*1000.0*{totLumi}*{xsec}/{weight_sum}".format(totLumi = TotalLumi, xsec = xsec1, weight_sum = event_weight_sum1)
    #weight2 = totalweight*"sample_weight*1000.0*{totLumi}*{xsec}/{weight_sum}".format(totLumi = TotalLumi, xsec = xsec2, weight_sum = event_weight_sum2)
    weight1 = "1";weight2="1"

    
    legs = ["Louvain","Mine"]
    for ch in channelcuts:
        chcut = channelcuts[ch]["cut"]
        cut1 = "("+chcut + " && " + llcut+")*"+ weight1
        cut2 = "("+chcut + " && " + llcut+")*"+ weight2
        #print "cut1 ",cut1, " cut2 ",cut2
        name = todraw
        if "*" in name:
            name = todraw.split("*")[0]
        plotname = os.path.join(variablesdir, name+"_"+ch)
        plotshapes(file1, treename1, todraw, cut1, file2, treename2, todraw, cut2, xbins, xtitle, legs, plotname)



def compareTH1(filename, histnamelist, legends, xbins, xtitle, plotname):
    colors = [800-4, 593, 820, 616-7, 432, 418]
    
    b1 = ROOT.TH1F("b1","b1",xbins[0],xbins[1],xbins[2])
    #b1.SetTitle("h2#rightarrow hh#rightarrow BBWW, B3"+" "*12 + "CMS Simulation Preliminary")
    b1.GetYaxis().SetTitle("Events")
    b1.GetXaxis().SetTitle("%s"%xtitle)
    b1.SetStats(0)
    histlist = []
    tf = ROOT.TFile(filename,"READ")
    hs = ROOT.THStack("Shapecomparison", "  ")
    legend = ROOT.TLegend(0.13,0.65,0.4,0.65+4*.04); 
    legend.SetFillStyle(0);
    legend.SetTextSize(0.04); legend.SetTextFont(42); legend.SetBorderSize(0)

    for i,histname in enumerate(histnamelist):
        hist = tf.Get(histname)
        ROOT.SetOwnership(hist, False)
        #hist.Print("all")
        hist.SetLineColor(colors[i])
        hist.SetLineWidth(2)
        #hs.Add(hist)
        histlist.append(hist)


    #for i,hist in enumerate(histlist):
        legend.AddEntry(hist, legends[i],"L")
        hs.Add(hist)
        #hist.Draw("sameL")

    c1 = ROOT.TCanvas("c1","c1",800,600)
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()
    #c1.Clear()
    #b1.Draw()
    hs.Draw("nostackL")
    hs.GetHistogram().GetXaxis().SetTitle(xtitle)
    hs.GetHistogram().GetYaxis().SetTitle("Events")
    tex0 = ROOT.TLatex(0.1,0.92, " #scale[1.4]{#font[61]{CMS}} Internal"+"  "*10+"35.87 fb^{-1} (13 TeV),2016")
    tex0.SetNDC(); tex0.SetTextSize(.045); tex0.SetTextFont(42)
    legend.Draw("same")
    tex0.Draw("same")


    c1.SaveAs(plotname+".pdf")
    c1.SaveAs(plotname+".C")

def compareTH1_v2(filenames, histnamelist, legends, filesuffix, xbins, xtitle, plotname):
    colors = [2, 593, 820, 616-7, 432, 418]
    style  = [1, 2]
    
    b1 = ROOT.TH1F("b1","b1",xbins[0],xbins[1],xbins[2])
    #b1.SetTitle("h2#rightarrow hh#rightarrow BBWW, B3"+" "*12 + "CMS Simulation Preliminary")
    b1.GetYaxis().SetTitle("Events")
    b1.GetXaxis().SetTitle("%s"%xtitle)
    b1.SetStats(0)
    #b1.SetMaximum(30000)
    #b1.Draw()
    histlist = []
    hs = ROOT.THStack("Shapecomparison", "  ")
    legend = ROOT.TLegend(0.2,0.65,0.4,0.65+4*.04); 
    legend.SetFillStyle(0);
    legend.SetTextSize(0.04); legend.SetTextFont(42); legend.SetBorderSize(0)

    tflist = []
    soversqrtb = []
    hs_soverb = ROOT.THStack("Shapecomparison_soverb", "  ")
    
    for j, filename in enumerate(filenames):
        tflist.append( ROOT.TFile.Open(filename))
        signal = None;  background = None;
        for i,histname in enumerate(histnamelist):
            hist = tflist[j].Get(histname)
            #ROOT.SetOwnership(hist, False)
            #hist.Print("all")
            #hist.SetLineColor(colors[i+j*len(histnamelist)])
            hist.SetLineColor(colors[j])
            hist.SetLineStyle(style[i])
            hist.SetLineWidth(2)
            hist.SetFillStyle(0)
            #hs.Add(hist)
            histlist.append(hist)
            #hist.Print("all")


        #for i,hist in enumerate(histlist):
            legend.AddEntry(hist, filesuffix[j]+" "+legends[i],"L")
            hs.Add(hist.Clone())
            #hist.Draw("sameL")
        signal = histlist[-1]; background = histlist [-2];
        soversqrtb.append(ROOT.TH1F(filesuffix[j],filesuffix[j],xbins[0],xbins[1],xbins[2]))
        for k in range(xbins[0]):
            s = signal.GetBinContent(k+1)
            b = background.GetBinContent(k+1)
            if b>0:
                soversqrtb[-1].SetBinContent(k+1, s/sqrt(b))
            else:
                soversqrtb[-1].SetBinContent(k+1, 100)
        soversqrtb[-1].SetLineColor(colors[j])
        soversqrtb[-1].SetLineWidth(2)
        hs_soverb.Add(soversqrtb[-1])

        

    c1 = ROOT.TCanvas("c1","c1",800,600)

    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(.0)
    pad1.Draw()
    pad1.cd()
    pad1.SetGridx()
    pad1.SetGridy()
    pad1.SetTickx()
    pad1.SetTicky()
    #pad1.SetLogy()
    #for hist in  histlist:
    #    hist.Print()
    #c1.Clear()
    #b1.Draw()
    hs.Draw("nostack,e")
    hs.GetHistogram().GetXaxis().SetTitle(xtitle)
    hs.GetHistogram().GetYaxis().SetTitle("Events")
    tex0 = ROOT.TLatex(0.1,0.92, " #scale[1.4]{#font[61]{CMS}} Internal"+"  "*10+"35.87 fb^{-1} (13 TeV),2016")
    tex0.SetNDC(); tex0.SetTextSize(.045); tex0.SetTextFont(42)
    legend.Draw("same")
    tex0.Draw("same")


    c1.cd()
    c1.Update()


    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, .29)
    pad2.SetTopMargin(0.)
    pad2.SetBottomMargin(.35)
    pad2.SetTicks(1,1)#draw x, y axis on both side (left right for y, and bottom up for x)
    pad2.SetGridy()
    pad2.Draw()
    pad2.cd()

    hs_soverb.Draw("nostackhist")
    hs_soverb.GetHistogram().GetXaxis().SetTitle(xtitle)
    hs_soverb.GetHistogram().GetYaxis().SetTitle("S/#sqrt{B}")
    hs_soverb.GetHistogram().GetXaxis().SetTitleSize(20)
    hs_soverb.GetHistogram().GetXaxis().SetTitleFont(43)
    hs_soverb.GetHistogram().GetXaxis().SetTitleOffset(3.0)
    hs_soverb.GetHistogram().GetXaxis().SetLabelSize(15)
    hs_soverb.GetHistogram().GetXaxis().SetLabelFont(43)#Absolute font size in pixel (precision 3)
    hs_soverb.GetHistogram().GetYaxis().SetNdivisions(505)
    hs_soverb.GetHistogram().GetYaxis().CenterTitle()
    hs_soverb.GetHistogram().GetYaxis().SetTitleSize(20)
    hs_soverb.GetHistogram().GetYaxis().SetTitleFont(43)
    hs_soverb.GetHistogram().GetYaxis().SetTitleOffset(.9)
    hs_soverb.GetHistogram().GetYaxis().SetLabelSize(15)
    hs_soverb.GetHistogram().GetYaxis().SetLabelFont(43)#Absolute font size in pixel (precision 3)



    c1.SaveAs(plotname+".pdf")
    c1.SaveAs(plotname+".C")



def compareKinematic_all(file1, file2, variablesdir):
    #jet1_pt = "jet1_pt*(jet1_pt>jet2_pt)+jet2_pt*(jet2_pt>jet1_pt)"
    #jet2_pt = "jet2_pt*(jet1_pt>jet2_pt)+jet1_pt*(jet2_pt>jet1_pt)"
    #jet1_eta = "jet1_eta*(jet1_pt>jet2_pt)+jet2_eta*(jet2_pt>jet1_pt)"
    #jet2_eta = "jet2_eta*(jet1_pt>jet2_pt)+jet1_eta*(jet2_pt>jet1_pt)"
    #compareKinematic(file1,file2, 'lep1_pt', [60, 10.0, 200], "lep1 p_{T}", variablesdir)
    #compareKinematic(file1,file2, 'lep2_pt', [60, 10.0, 200], "lep2 p_{T}", variablesdir)
    #compareKinematic(file1,file2, 'jet1_pt', [70, 20.0, 300], "jet1 p_{T}", variablesdir)
    #compareKinematic(file1,file2, 'jet2_pt', [70, 20.0, 300], "jet2 p_{T}", variablesdir)
    #compareKinematic(file1,file2, 'lep1_eta', [60, -2.4, 2.4], "lep1 #eta", variablesdir)
    #compareKinematic(file1,file2, 'lep2_eta', [60, -2.4, 2.4], "lep2 #eta", variablesdir)
    #compareKinematic(file1,file2, 'jet1_eta', [70, -2.5, 2.5], "jet1 #eta", variablesdir)
    #compareKinematic(file1,file2, 'jet2_eta', [70, -2.5, 2.5], "jet2 #eta", variablesdir)
    #compareKinematic(file1,file2, 'met_pt', [50, 0.0, 500.0],"MET p_{T}", variablesdir)
    #compareKinematic(file1,file2, 'met_phi', [60, -3.2, 3.20],"MET #phi", variablesdir)
    compareKinematic(file1,file2, 'll_M', [50, 12.0, 100.0], "M_{ll}", variablesdir)
    compareKinematic(file1,file2, 'll_pt', [50, 0.0, 450.0], "Dilepton p_{T}", variablesdir)
    compareKinematic(file1,file2, 'jj_pt', [50, 0.0, 450.0], "Dijet p_{T}", variablesdir)
    #compareKinematic(file1,file2, 'll_M', [50, 12.0, 76.0], "M_{ll}", variablesdir)
    #compareKinematic(file1,file2, 'll_DR_l_l', [50, .0, 6.0], "#DeltaR_{ll}", variablesdir)
    #compareKinematic(file1,file2, 'jj_M', [50, 0.0, 400.0], "M_{jj}",variablesdir)
    #compareKinematic(file1,file2, 'jj_DR_j_j', [50, .0, 6.0], "#DeltaR_{jj}",variablesdir)
    #compareKinematic(file1,file2, 'llmetjj_DPhi_ll_jj', [24, .0, 3.1415926],"#Delta#phi(ll,jj)", variablesdir)
    #compareKinematic(file1,file2, 'll_pt', [50, 0.0, 450.0], "Dilepton p_{T}", variablesdir)
    #compareKinematic(file1,file2, 'jj_pt', [50, 0.0, 450.0], "Dijet p_{T}", variablesdir)
    #compareKinematic(file1,file2, 'llmetjj_minDR_l_j', [50, .0, 5.0], "#DeltaR_{l,j}", variablesdir)
    #compareKinematic(file1,file2, 'llmetjj_MTformula', [50, 0.0, 500.0],"MT", variablesdir)
    #compareKinematic(file1,file2, 'met_pt', [50, 0.0, 500.0],"MET", variablesdir)



def compareSignalAndBackground(masslist, filepath, DNNtype, outdir):
    xtitle = "DNN output"
    xbins = [75, 0.0, 3.0]
    ##Hhh_FinalBGYield_xsec1pb_NN_nnout_MTonly_nnstep0p04_nncut0p72_SignalM900
    allChannels = ["MuMu", "ElEl", "MuEl"]
    for i, mass in enumerate(masslist):
        plotname = "CompareFinalshapes_%s"%(DNNtype)
        histnamelist = []; legends = []
        filename = "Hhh_FinalBGYield_xsec1pb_NN_nnout_%s_nnstep0p04_nncut0p0_SignalM%d.root"%(DNNtype, mass)
        filename = os.path.join(filepath, filename)
        for channel in allChannels:
            thispath = os.path.join(filepath, filename)
            histnamelist.append("bg_all_%s_M%d"%(channel, mass))
            histnamelist.append("RadionM%d_%s"%(mass, channel))
            legends.append("%s_Background_M%d"%(channel, mass))
            legends.append(" %s_Signal_M%d"%(channel, mass))
        plotname = plotname+"_Mass%d"%mass
        plotname = os.path.join(outdir, plotname)
        print  " histnames ",histnamelist," plotname ", plotname


        compareTH1(filename, histnamelist, legends, xbins, xtitle, plotname)

    
def compareSignalAndBackground_v2(masslist, filepath, filesuffix, DNNtype, outdir):
    xtitle = "DNN output"
    #xbins = [25, 0.0, 1.0]
    xbins = [75, 0.0, 3.0]
    ##Hhh_FinalBGYield_xsec1pb_NN_nnout_MTonly_nnstep0p04_nncut0p72_SignalM900
    allChannels = ["MuMu", "ElEl", "MuEl"]
    for i, mass in enumerate(masslist):
        filename1 = "Hhh_FinalBGYield_xsec1pb_NN_nnout_%s_nnstep0p04_nncut0p0_SignalM%d.root"%(DNNtype, mass)
        filename1 = os.path.join(filepath[0], filename1)
        filename2 = "Hhh_FinalBGYield_xsec1pb_NN_nnout_%s_nnstep0p04_nncut0p0_SignalM%d.root"%(DNNtype, mass)
        filename2 = os.path.join(filepath[1], filename2)
        filenames = [filename1, filename2]
        #filenames = [filename1]
        for channel in allChannels:
            histnamelist = []; legends = []
            
            #thispath = os.path.join(filepath, filename)
            #histnamelist.append("bg_all_%s_M%d"%(channel, mass))
            histnamelist.append("sT_%s"%(channel))
            histnamelist.append("RadionM%d_%s"%(mass, channel))
            #legends.append("%s_Background_M%d"%(channel, mass))
            legends.append("%s_sT"%(channel))
            legends.append(" %s_Signal_M%d"%(channel, mass))
            plotname = "CompareFinalshapes_%s"%(DNNtype)+"_"+channel+"_Mass%d"%mass
            plotname = os.path.join(outdir, plotname)
            print  " histnames ",histnamelist," plotname ", plotname
            compareTH1_v2(filenames, histnamelist, legends, filesuffix,xbins, xtitle, plotname)

    



outputdir = "Louvain_Mine_weight_comparison_noPUweight_newTT/"
os.system("mkdir -p "+outputdir)
file_L = "/Users/taohuang/Documents/DiHiggs/20170530/HHNTuples/TTTo2L2Nu_13TeV-powheg_Summer16MiniAODv2_v5.0.1+80X_HHAnalysis_2017-03-01.v0_histos.root"

file_Tao = "TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8_Friend.root"
#compareLouvainandMine(file_L, file_Tao, outputdir)

variablesdir = "Louvain_Mine_kinematics_comparison_noweight_newTT/"

#os.system("mkdir -p "+variablesdir)
#compareKinematic_all(file_L, file_Tao, variablesdir)
masslist = [400, 900]
masslist = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900]
outdir = "CompareSignalandBackgroundShapes_TT_StatisticalError/"
filepath = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy1D/"
os.system("mkdir -p "+outdir)
#compareSignalAndBackground(masslist, filepath, "MTonly", outdir)

filesuffix = ["Merged","SignalRegionOnly"]
filesuffix = ["Old","New"]
filepaths = ["HHbbWW_20200623_NNoutput_NoMjjCR_NNcutstudy1D/", "HHbbWW_20200623_NNoutput_NoMjjCRMjjcut_NNcutstudy1D/"]
filepaths = ["HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy1D/","HHbbWW_20200807_NNoutput_MjjCR_NNcutstudy1D/"]
compareSignalAndBackground_v2(masslist, filepaths, filesuffix, "MTandMT2_MJJ", outdir)
