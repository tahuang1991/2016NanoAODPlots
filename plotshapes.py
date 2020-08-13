import os
import numpy as np
import ROOT



#colors = [623, 593, 820, 797, 432, 418]
colors = [800-4, 593, 820, 616-7, 432, 418]
bghistnamedict = {"TT":"ttbar","ttV":"ttV", "singleTop":"SingleTop","SMH":"SMHiggs","VV":"VV","Drell-Yan":"dy_mc"}
def plotshapes(filedir, mass, bgnames, plotname):
    
    xtitle = "NN output, Mjj bin"
    suffix = "CMSHIG17006"
    channels = ["ElEl","MuEl","MuMu"]
    for ich, channel in enumerate(channels):
        allhists = []
        rfile = ROOT.TFile(filedir+"GGToX0ToHHTo2B2L2Nu_M%d_%s_shapes.root"%(mass, channel), "READ")
        hs = ROOT.THStack("allbg_"+channel+"_"+suffix, "  ")
        legend = ROOT.TLegend(0.7,0.65,0.88,0.65+len(bgnames)*.045); 
	legend.SetTextSize(0.04); legend.SetTextFont(42); legend.SetBorderSize(0)
        legend2 = ROOT.TLegend(0.5,0.75,0.68,0.88); 
	legend2.SetTextSize(0.04); legend2.SetTextFont(42); legend2.SetBorderSize(0)
        #legend.SetHeader("DNN training: kinematics+")
	BGSum = 0.0
        maxbgbin = 0.0
        hist_data = rfile.Get("mjj_vs_NN_M%d/data_obs"%mass)
	for i, key in enumerate(bgnames):
            histname  =  bghistnamedict[key]
            if channel != "MuEl" and key == "Drell-Yan":
                histname = "data_nobtag_to_btagM"
            hist = rfile.Get("mjj_vs_NN_M%d/"%mass + histname) 
            hist.SetStats(0)
            if key == "Drell-Yan" and channel != "MuEl":
                hist1 = rfile.Get("mjj_vs_NN_M%d/"%mass + "ttbar_nobtag_to_btagM") 
                hist2 = rfile.Get("mjj_vs_NN_M%d/"%mass + "SingleTop_nobtag_to_btagM") 
                hist.Add(hist1, -1)
                hist.Add(hist2, -1)
	    print "mass ",mass, " channel ", channel," bg ",key," rate ",hist.Integral()," histname ",histname
	    hs.Add(hist)
	    hist.SetFillColor(colors[i])
	    allhists.append(hist)
	    BGSum = BGSum + hist.Integral()
        print "allhists  ", allhists
        hist_bgall = allhists[0].Clone()
	for i in range(len(allhists)):
            if i>0:
                hist_bgall.Add(allhists[i])
	    legend.AddEntry(allhists[len(allhists)-i-1], bgnames[len(allhists)-i-1], "f")
	
        #hist_s = ROOT.TH1F("signal_"+channel+"_%s"%(suffix), "signal_"+channel+"_%s"%( suffix), len(xbins)-1, xbins)
        hist_s = rfile.Get("mjj_vs_NN_M%d/ggX0HH%d"%(mass, mass))
        hist_s.Scale(5000.0)#1fb->5pb
        maxbgbin = hist_data.GetBinContent(hist_data.GetMaximumBin())
        maxsignalbin = hist_s.GetBinContent(hist_s.GetMaximumBin())
        if maxsignalbin>maxbgbin:
            hs.SetMaximum(maxsignalbin*1.4)
        else:
            hs.SetMaximum(maxbgbin*1.3)

	print "mass ",mass, " channel ", channel," rate:signal ",hist_s.Integral()," BG(all) ",BGSum," data ",hist_data.Integral()

	hist_s.SetLineColor(ROOT.kRed)
	hist_s.SetLineWidth(2)
        hist_data.SetMarkerStyle(20)
        hist_data.SetMarkerColor(ROOT.kBlack)
        legend2.AddEntry(hist_data, "Data","p")
        legend2.AddEntry(hist_data," ","")
        signaltext = "#splitline{Signal(5pb)}{M=%d GeV}"%(mass)
        legend2.AddEntry(hist_s, signaltext,"l")
        #legend2.AddEntry(hist_s, "Signal, 5pb","l")




        c1 = ROOT.TCanvas("c", "canvas", 800, 800)
        c1.Clear()
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(.0)
        pad1.Draw()
        pad1.cd()
	hs.Draw("hist")
	hist_s.Draw("samehist")
        hist_data.Draw("samepe")
	
        #tex1 = ROOT.TLatex(0.17,0.8, channelcuts[channel]["latex"]+" channel, "+nnout)
	tex1 = ROOT.TLatex(0.17,0.8, channel+ " channel ")
	tex1.SetNDC(); tex1.SetTextSize(.045)
	tex2 = ROOT.TLatex(0.15,0.56, "M_{jj}<75 GeV "+"  "*8+" 75 GeV <=M_{jj}<140 GeV"+ "  "*8+" M_{jj} >= 140 GeV")
	tex2.SetNDC(); tex2.SetTextSize(.04)
	tex1.Draw("same")
        tex2.Draw("same")
	#hs.GetHistogram().SetTitle(" #scale[1.4]{#font[61]{CMS}} Simulation Preliminary"+"  "*38+"35.87 fb^{-1} (13 TeV),2016")
        #hs.GetHistogram().SetTitle("")
	#hs.GetHistogram().SetTitleSize(.04)
	#hs.GetHistogram().SetTitleOffset(1.2)
	tex0 = ROOT.TLatex(0.1,0.92, " #scale[1.4]{#font[61]{CMS}} "+"  "*24+"35.87 fb^{-1} (13 TeV),2016")
	tex0.SetNDC(); tex0.SetTextSize(.05); tex0.SetTextFont(42)
	tex0.Draw("same")
        #hs.GetHistogram().GetXaxis().SetTitle("DNN output, M_{jj} bins")
        #hs.GetHistogram().GetXaxis().SetTitle(xtitle)
        hs.GetHistogram().GetYaxis().SetTitle("Events")
        legend.Draw("same")
        legend2.Draw("same")
	#tex.Draw("same")
        #add data vs MC
        c1.cd()
        c1.Update()
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, .29)
        pad2.SetTopMargin(0.)
        pad2.SetBottomMargin(.35)
        pad2.SetTicks(1,1)#draw x, y axis on both side (left right for y, and bottom up for x)
        pad2.SetGridy()
        pad2.Draw()
        pad2.cd()


        hist_ratio = hist_data.Clone()
        hist_ratio.Divide(hist_bgall)

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
        hist_ratio.GetYaxis().SetTitle("Data/MC")
        hist_ratio.GetYaxis().SetNdivisions(505)
        hist_ratio.GetYaxis().CenterTitle()
        hist_ratio.GetYaxis().SetTitleSize(20)
        hist_ratio.GetYaxis().SetTitleFont(43)
        hist_ratio.GetYaxis().SetTitleOffset(.9)
        hist_ratio.GetYaxis().SetLabelSize(15)
        hist_ratio.GetYaxis().SetLabelFont(43)#Absolute font size in pixel (precision 3)


        c1.SaveAs(plotname+"_Radion_"+channel+"_"+suffix+".C")
        #c1.SaveAs(plotname+"_Radion_"+channel+".png")
        c1.SaveAs(plotname+"_Radion_"+channel+"_"+suffix+".pdf")


bghistnamedict = {"TT":"ttbar","ttV":"ttV", "sT":"SingleTop","SMH":"SMHiggs","VV":"VV","DY":"dy_mc"}
#extraweight is used for compensating the "HLT" safe ID cut, for Louvain ntuples
channelcuts = {"MuMu":{"cut":"isMuMu","extraweight": 24827.9/23891.14, "latex":"#mu#mu"},
		"MuEl":{"cut":"(isMuEl || isElMu)", "extraweight": 25786.8/30211.9, "latex":"e#mu"},
		"ElEl":{"cut":"isElEl","extraweight":8487.6/12190.1, "latex":"ee"},
		}
def CompareShapes(filedir, filesignal, backgrounddict, todraw, cut, xbins, mass, bgnames, xtitle, suffix, plotname):
    
    #if len(xbins) == 3:
    #    nbins = xbins[0]; xmin = xbins[1]; xmax =  xbins[2]
    #    xbins = []
    #    binwidth = (xmax-xmin)/nbins
    #    for i in range(0, nbins+1):
    #        xbins.append(xmin + i*binwidth)
    #    xbins = np.asarray(xbins)

    #xtitle = "NN output, Mjj bin"
    #suffix = "CMSHIG17006"
    channels = ["ElEl","MuEl","MuMu"]

    chlist = {}
    treename = "evtreeHME_nn"
    #treename = "t"
    for key in backgrounddict:
	ch = ROOT.TChain(treename)
	for f in backgrounddict[key]:
	    ch.Add(f)
	chlist[key] = ch
    ch_s =  ROOT.TChain(treename); ch_s.Add(filesignal)
    xsec_signal = 5.0#pb
    weight = "(35.87*1000*cross_section/event_weight_sum)*final_total_weight"
    weight_s = "(35.87*1000*{xsec}/event_weight_sum)*final_total_weight".format(xsec = xsec_signal)

    for ich, channel in enumerate(channels):
        allhists = []
        rfile = ROOT.TFile(filedir+"GGToX0ToHHTo2B2L2Nu_M%d_%s_shapes.root"%(mass, channel), "READ")
        legend = ROOT.TLegend(0.7,0.7,0.88,0.7+4*.045); 
        legend.SetHeader("Official")
	legend.SetTextSize(0.04); legend.SetTextFont(42); #legend.SetBorderSize(0)
        legend2 = ROOT.TLegend(0.45,0.7+0.045,0.62,0.7+4*.045); 
	legend2.SetTextSize(0.04); legend2.SetTextFont(42); #legend2.SetBorderSize(0)
        legend2.SetHeader("Mine")
        
        #legend.SetHeader("DNN training: kinematics+")
	BGSum = 0.0
        maxbgbin = 0.0
        #### hist from official ntuples
        hist_data = rfile.Get("mjj_vs_NN_M%d/data_obs"%mass)
    
	for i, key in enumerate(bgnames):
            histname  =  bghistnamedict[key]
            #if channel != "MuEl" and key == "Drell-Yan":
            if channel != "MuEl" and key == "DY":
                histname = "data_nobtag_to_btagM"
            hist = rfile.Get("mjj_vs_NN_M%d/"%mass + histname) 
            hist.SetStats(0)
            #if key == "Drell-Yan" and channel != "MuEl":
            if key == "DY" and channel != "MuEl":
                hist1 = rfile.Get("mjj_vs_NN_M%d/"%mass + "ttbar_nobtag_to_btagM") 
                hist2 = rfile.Get("mjj_vs_NN_M%d/"%mass + "SingleTop_nobtag_to_btagM") 
                hist.Add(hist1, -1)
                hist.Add(hist2, -1)
	    print "mass ",mass, " channel ", channel," bg ",key," rate ",hist.Integral()," histname ",histname
            #hist.SetFillColor(colors[i])
	    allhists.append(hist)
	    BGSum = BGSum + hist.Integral()
        print "allhists  ", allhists
        hist_bgall = allhists[0].Clone()
	for i in range(len(allhists)):
            if i>0:
                hist_bgall.Add(allhists[i])
            #legend.AddEntry(allhists[len(allhists)-i-1], bgnames[len(allhists)-i-1], "f")
	

        #hist_s = ROOT.TH1F("signal_"+channel+"_%s"%(suffix), "signal_"+channel+"_%s"%( suffix), len(xbins)-1, xbins)
        hist_s = rfile.Get("mjj_vs_NN_M%d/ggX0HH%d"%(mass, mass))
        hist_s.Scale(xsec_signal*1000.0)#1fb->5pb
        maxbgbin = hist_data.GetBinContent(hist_data.GetMaximumBin())
        maxsignalbin = hist_s.GetBinContent(hist_s.GetMaximumBin())

        
	hist_s.SetLineColor(ROOT.kRed)
	hist_s.SetLineWidth(2)
        hist_data.SetMarkerStyle(20)
        hist_data.SetMarkerColor(ROOT.kBlack)
        hist_bgall.SetLineColor(ROOT.kBlue)
	hist_bgall.SetLineWidth(2)
        
	print "mass ",mass, " channel ", channel," from official, rate:signal ",hist_s.Integral()," BG(all) ",BGSum," data ",hist_data.Integral()

        legend.AddEntry(hist_data, "Data","p")
        legend.AddEntry(hist_bgall, "Backgrounds","l")
        #signaltext = "#splitline{Signal(5pb)}{M=%d GeV}"%(mass)
        signaltext = "Signal(5pb)"
        legend.AddEntry(hist_s, signaltext,"l")

        ### from evtree
        allhists_v2 = []
        #hist_data_fake = ROOT.TH1F("data_obs_"+channel+"_M%d_%s"%(mass, suffix), "data_obs_"+channel+"_M%d_%s"%(mass, suffix), len(xbins)-1, xbins)
	hist_data_fake = ROOT.TH1F("data_obs_"+channel+"_%s"%(suffix), "data_obs_"+channel+"_%s"%( suffix), xbins[0],xbins[1], xbins[2])
	for i, key in enumerate(bgnames):
	    hist = ROOT.TH1F(key+"_"+channel+"_%s"%(suffix), key+"_"+channel+"_%s"%(suffix),  xbins[0],xbins[1], xbins[2])
            finalcut = "("+ cut + " && "+ channelcuts[channel]["cut"] +")*"+weight+"*%f"%(channelcuts[channel]["extraweight"])
            #finalcut = "("+ cut + " && "+ channelcuts[channel]["cut"] +")*"+weight
            print "todraw ",todraw," finalcut ",finalcut
            hist.SetStats(0)
	    chlist[key].Draw(todraw + ">> " + hist.GetName(), finalcut)
	    hist.SetFillColor(colors[i])
	    allhists_v2.append(hist)
            #print "mass ",mass, " channel ", channel," bg ",key," rate ",hist.Integral()
            #if key == "TT":
            #    maxbgbin = hist.GetBinContent(hist.GetMaximumBin())
	    hist_data_fake.Add(hist)
	#for hist,name in zip(reversed(allhists), reversed(bgnames)):
        #for i in range(len(allhists_v2)):
            #legend.AddEntry(allhists_v2[len(allhists)-i-1], bgnames[len(allhists)-i-1], "f")
	
        BGSum_v2 = hist_data_fake.Integral()
        #hist_s_v2 = ROOT.TH1F("signal_"+channel+"_%s"%(suffix), "signal_"+channel+"_%s"%( suffix), len(xbins)-1, xbins)
        hist_s_v2 = ROOT.TH1F("signal_"+channel+"_%s"%(suffix), "signal_"+channel+"_%s"%( suffix),  xbins[0],xbins[1], xbins[2])
        finalcut_s  = "("+ cut + " && "+ channelcuts[channel]["cut"] +")*"+weight_s+"*%f"%(channelcuts[channel]["extraweight"])
        #finalcut_s  = "("+ cut + " && "+ channelcuts[channel]["cut"] +")*"+weight_s
	ch_s.Draw(todraw + ">> " + hist_s_v2.GetName(), finalcut_s)

	print "mass ",mass, " channel ", channel," from evtree, rate:signal ",hist_s_v2.Integral()," BG(TT,DY,sT) ",BGSum_v2," data(fake) ",hist_data_fake.Integral()

	hist_s_v2.SetLineColor(ROOT.kRed)
	hist_s_v2.SetLineWidth(2)
	hist_s_v2.SetLineStyle(2)
	hist_data_fake.SetLineColor(ROOT.kBlue)
	hist_data_fake.SetLineWidth(2)
	hist_data_fake.SetLineStyle(2)
        maxbgbin_v2 = hist_data_fake.GetBinContent(hist_data_fake.GetMaximumBin())
        maxsignalbin_v2 = hist_s_v2.GetBinContent(hist_s_v2.GetMaximumBin())


        legend2.AddEntry(hist_data_fake, "Backgrounds","l")
        #signaltext = "#splitline{Signal(5pb)}{M=%d GeV}"%(mass)
        signaltext = "Signal(5pb)"
        legend2.AddEntry(hist_s_v2, signaltext,"l")
        #legend2.AddEntry(hist_s, "Signal, 5pb","l")




        c1 = ROOT.TCanvas("c", "canvas", 800, 800)
        c1.Clear()
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(.0)
        pad1.Draw()
        pad1.cd()
        if maxbgbin < maxbgbin_v2:
            hist_bgall.SetMaximum(maxbgbin_v2 * 1.2)
        else:
            hist_bgall.SetMaximum(maxbgbin*1.2)
        hist_bgall.SetTitle(" ")
        hist_bgall.Draw("histe2")
	hist_s.Draw("samehist")
        hist_data.Draw("sameep")

        hist_data_fake.Draw("histsame")
        hist_s_v2.Draw("histsame")
	
        
        #tex1 = ROOT.TLatex(0.17,0.8, channelcuts[channel]["latex"]+" channel, "+nnout)
	tex1 = ROOT.TLatex(0.17,0.8, channel+ " channel, "+suffix)
	tex1.SetNDC(); tex1.SetTextSize(.045)
	tex2 = ROOT.TLatex(0.18,0.56, "M_{jj}<75 GeV "+"  "*10+" 75 GeV <=M_{jj}<140 GeV"+ "  "*12+" M_{jj} >= 140 GeV")
	tex2.SetNDC(); tex2.SetTextSize(.04)
	tex1.Draw("same")
        tex2.Draw("same")
	tex0 = ROOT.TLatex(0.1,0.92, " #scale[1.4]{#font[61]{CMS}} "+"  "*30+"35.87 fb^{-1} (13 TeV),2016")
	tex0.SetNDC(); tex0.SetTextSize(.05); tex0.SetTextFont(42)
	tex0.Draw("same")
        hist_bgall.GetXaxis().SetTitle(xtitle)
        legend.Draw("same")
        legend2.Draw("same")
	#tex.Draw("same")
        #add data vs MC
        c1.cd()
        c1.Update()
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, .29)
        pad2.SetTopMargin(0.)
        pad2.SetBottomMargin(.35)
        pad2.SetTicks(1,1)#draw x, y axis on both side (left right for y, and bottom up for x)
        pad2.SetGridy()
        pad2.Draw()
        pad2.cd()


        hist_ratio = hist_data_fake.Clone()
        hist_ratio.Divide(hist_bgall)
        hist_ratio.SetMarkerColor(ROOT.kBlue)


        hist_ratio_s = hist_s_v2.Clone()
        hist_ratio_s.Divide(hist_s)
        hist_ratio_s.SetMarkerColor(ROOT.kRed)
        print "GetMean(), ratio bg ",hist_ratio.GetMean()," signal ",hist_ratio_s.GetMean()
        ratio_entries_bg = hist_data_fake.Integral()*1.0/hist_bgall.Integral()
        ratio_entries_s = hist_s_v2.Integral()*1.0/hist_s.Integral()
        print "ratio bg ",ratio_entries_bg," signal ",ratio_entries_s


        hist_ratio.Draw("ep")
        hist_ratio.SetStats(0)
        hist_ratio.SetTitle("")
        hist_ratio_s.Draw("epsame")
        hist_ratio_s.SetStats(0)
        hist_ratio_s.SetTitle("")
        maxratiobin = hist_ratio.GetBinContent(hist_ratio.GetMaximumBin())
        minratiobin = hist_ratio.GetBinContent(hist_ratio.GetMinimumBin())
        #maxratiobin_s = hist_ratio.GetBinContent(hist_ratio_s.GetMaximumBin())
        #minratiobin_s = hist_ratio.GetBinContent(hist_ratio_s.GetMinimumBin())
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
        hist_ratio.GetYaxis().SetTitle("Mine/Official")
        hist_ratio.GetYaxis().SetNdivisions(505)
        hist_ratio.GetYaxis().CenterTitle()
        hist_ratio.GetYaxis().SetTitleSize(20)
        hist_ratio.GetYaxis().SetTitleFont(43)
        hist_ratio.GetYaxis().SetTitleOffset(.9)
        hist_ratio.GetYaxis().SetLabelSize(15)
        hist_ratio.GetYaxis().SetLabelFont(43)#Absolute font size in pixel (precision 3)
        tex = ROOT.TLatex(0.15,0.45, "Overall ratio, signal: %.2f, BG: %.2f"%(ratio_entries_s, ratio_entries_bg))
        tex.SetNDC()
        tex.SetTextSize(.07)
        tex.Draw("same")


        c1.SaveAs(plotname+"_Radion_"+channel+"_"+suffix+".C")
        #c1.SaveAs(plotname+"_Radion_"+channel+".png")
        c1.SaveAs(plotname+"_Radion_"+channel+"_"+suffix+".pdf")


def CompareLouvainShapehist(masspoints, nnout, output_folder, outdir):
    bgnames = ["sT","DY","TT"]
    #bgnames = ["TT"]
    plotname = os.path.join(outdir, "Louvainshapes_xsec5pb_NN_%s"%nnout)
    ###create tfile
    nnbins = []
    for x in range(0, 75-9+1):
        nnbins.append((x+3)*1.0/25.0)
    print "nnbins ",nnbins
    #nnbins_x = np.asarray(nnbins)
    nnbins_x = [66, 3./25, 3-6.0/25]
    xtitle = "DNN output, M_{jj} bins"
    for mass in masspoints:
        filedir = "M%d.r7526/"%mass
        #cut = "({nnout}_M{mass}>3.0/25 && hme_h2mass_reco>=250)".format(nnout = nnout, mass = mass)
        cut = "({nnout}_M{mass}>3.0/25 )".format(nnout = nnout, mass = mass)
        todraw = "(({nnout}_M{mass})*(jj_M<75 && {nnout}_M{mass}>3.0/25)+(jj_M>=75 && jj_M<140 && {nnout}_M{mass}>3.0/25)*({nnout}_M{mass}+1-3.0/25)+(jj_M>=140 && {nnout}_M{mass}>3.0/25)*({nnout}_M{mass}+2-6.0/25))".format(nnout = nnout, mass=mass)
	file_s = os.path.join(output_folder, "radion_M%d_addNN.root"%(mass))
        backgrounddict = {}
	for bg in bgnames:
	    allfiles = os.listdir(output_folder)
	    backgrounddict[bg] = []
	    for f in allfiles:
		if f.startswith(bg):
		   backgrounddict[bg].append(os.path.join(output_folder, f))

        print "file_s ",file_s," file backgrounds ",backgrounddict, " official filedir ",filedir
        suffix = "M%d"%(mass)
        CompareShapes(filedir, file_s, backgrounddict, todraw, cut, nnbins_x, mass, bgnames, xtitle, suffix, plotname)
        



        
bgnames = ["TT","Drell-Yan","singleTop","ttV","SMH"]
mass = 400
outdir = "Shapeplots_Comparison/"
os.system("mkdir -p "+ outdir)
masspoints = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900]
for mass in masspoints:
    filedir = "M%d.r7526/"%mass
    plotname = os.path.join(outdir, "final_NN_Mjj_shape_M%d"%mass)
    #plotshapes(filedir, mass, bgnames, plotname)

output_folder = "/Users/taohuang/Documents/DiHiggs/20180126/20180205_20180202_10k_Louvain_ALL/2018-02-07_addNN/"
#masspoints = [400]
CompareLouvainShapehist(masspoints, "nnout_MTonly", output_folder, outdir)
