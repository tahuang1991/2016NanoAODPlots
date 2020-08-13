import os
from array import array
import datetime
import keras
from math import sqrt
import numpy.lib.recfunctions as recfunctions
import numpy as np
import ROOT

import sys
sys.argv.append( '-b' )
sys.argv.append( '-q' )

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

##local 
from common import *
import SampleTypeLUT


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.Reset()


#ROOT.gROOT.SetBatch(1)
#ROOT.gStyle.SetStatW(0.07)
#ROOT.gStyle.SetStatH(0.06)


#ROOT.gStyle.SetErrorX(0)
#ROOT.gStyle.SetErrorY(0)

#ROOT.gStyle.SetTitleStyle(0)
#ROOT.gStyle.SetTitleAlign(13) ## coord in top left
#ROOT.gStyle.SetTitleX(0.)
#ROOT.gStyle.SetTitleY(1.)
#ROOT.gStyle.SetTitleW(1)
#ROOT.gStyle.SetTitleH(0.058)
#ROOT.gStyle.SetTitleBorderSize(0)

#ROOT.gStyle.SetPadLeftMargin(0.126)
#ROOT.gStyle.SetPadRightMargin(0.10)
#ROOT.gStyle.SetPadTopMargin(0.06)
#ROOT.gStyle.SetPadBottomMargin(0.13)

#ROOT.gStyle.SetMarkerStyle(1)
### expected entries and new weight should be added if some subjobs are failed
expectedEntries = {'DYToLL1J': 5095747, 'radion_M600': 25970, 'sT_top': 27069, 'radion_M260': 9679, 'radion_M750': 29187, 'TT': 5487534, 'radion_M270': 10015, 'radion_M400': 18036, 'radion_M650': 27518, 'radion_M500': 22384, 'radion_M350': 15166, 'DYM10to50': 69311, 'radion_M450': 20295, 'radion_M800': 30349, 'radion_M900': 31338, 'radion_M300': 12251, 'DYToLL2J': 15265923, 'sT_antitop': 27460, 'DYToLL0J': 838304, 'radion_M550': 24578}

def plotmodel(modelLUT, model, plotname):
    modelfile = os.path.join(modelLUT[model]['workingdir'], "hh_resonant_trained_model.h5")
    thismodel = keras.models.load_model(modelfile)
    keras.utils.plot_model(thismodel, show_shapes = True, to_file=plotname)


def writeNNToTree(input_file, weight, modellist, modelLUT, sampletype, masslist, outFile):
    isSignal = (sampletype<13 or sampletype>=100)
    cut = "(91 - ll_M) > 15.0"
    FullInputs = ['jj_pt', 'll_pt', 'll_M', 'll_DR_l_l', 'jj_DR_j_j', 'llmetjj_DPhi_ll_jj', 'llmetjj_minDR_l_j', 'llmetjj_MTformula','mt2', 'jj_M','isSF','hme_h2mass_reco','isElEl','isElMu','isMuEl','isMuMu','met_pt','cross_section','event_weight_sum','sample_weight','event_reco_weight']
    #FullInputs = ['jj_pt', 'll_pt', 'll_M', 'll_DR_l_l', 'jj_DR_j_j', 'llmetjj_DPhi_ll_jj', 'llmetjj_minDR_l_j', 'llmetjj_MTformula', 'jj_M','isSF','isElEl','isElMu','isMuEl','isMuMu','met_pt','sample_weight','event_reco_weight']
    #dataset, weight = tree_to_numpy(f, FullInputs, weight_expression, cut, reweight_to_cross_section=False)
    file_handle = TFile.Open(input_file)

    #tree = file_handle.Get('t')
    tree = file_handle.Get('Friends')

    """
    ## not normalize to cross section
    cross_section = 1
    relative_weight = 1
    if not isSignal:
        print input_file
	relative_weight = getReweight_to_xsec(tree)
	print "reweight ",relative_weight

    if isinstance(weight, dict):
        # Keys are regular expression and values are actual weights. Find the key matching
        # the input filename
        found = False
        weight_expr = None
        if '__base__' in weight:
            weight_expr = weight['__base__']

        for k, v in weight.items():
            if k == '__base__':
                continue

            groups = re.search(k, input_file)
            if not groups:
                continue
            else:
                if found:
                    raise Exception("The input file is matched by more than one weight regular expression. %r" % input_file)

                found = True
                weight_expr = join_expression(weight_expr, v)

        if not weight_expr:
            raise Exception("Not weight expression found for input file %r" % weight_expr)

        weight = weight_expr

    """
    # Read the tree and convert it to a numpy structured array
    #a = tree2array(tree, branches= FullInputs + [weight], selection=cut)
    a = tree2array(tree, branches= FullInputs, selection=cut)
    # Rename the last column to 'weight'
    #print "print a type ", type(a)," shape ",a.shape, " len(a) ",len(a), " a.dtype.names ",a.dtype.names, " a[:10,:] ", a[:10]
    a.dtype= {"names": FullInputs ,"formats":[np.double]*(len(FullInputs)+1)}

  
    #print "print a type ", type(a)," shape ",a.shape, " a.dtype.names ",a.dtype.names, " a[:10,:] ", a[:10]
    
    
    ## add mass column
    """
    resonant_mass_probabilities = [[260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 800, 900], [ 0.0384008 ,  0.03985746,  0.04913783,  0.0611964 ,  0.07281586, 0.08038202,  0.09056594,  0.09963943,  0.10548293,  0.11182466, 0.12333857,  0.12735809]]

    masscol = np.array([0.0]*len(a))
    if isSignal:
	mass = SampleTypeLUT.getSignalMass(sampletype)
    	masscol  = np.array([mass]*len(a))
    else:
	#rs = np.random.RandomState(42)
	#indices = rs.choice(len(resonant_mass_probabilities[0]), len(a), p=resonant_mass_probabilities[1])
	#masscol = np.array(np.take(resonant_mass_probabilities[0], indices, axis=0), dtype='float')
	
	masscol = np.array([mass]*len(a))
    """
    allnnout = [] 
    nnoutname = []
    
    #sample_type = ROOT.vector('string')()
    for i, key in enumerate(modellist):
	modelfile = os.path.join(modelLUT[key]['workingdir'], "hh_resonant_trained_model.h5")
	inputs = modelLUT[key]['inputs'] 
        print "inputs ", inputs, " len ", len(inputs), " modelfile ",modelfile
	thismodel = keras.models.load_model(modelfile)
        for mass in masslist:	
	    thisdataset = np.array(a[inputs].tolist(), dtype=np.float32)
            print "thisdataset ",thisdataset.shape,"  ",thisdataset[:2]
	    masscol = [mass]*len(a)
	    thisdataset = np.c_[thisdataset, masscol]
	    #print "this model: ",key," inputs ",inputs, " datasetshape ",thisdataset.shape," type ",type(thisdataset)
	    nnout = thismodel.predict(thisdataset, batch_size=5000, verbose=1)
	    nnout = np.delete(nnout, 1, axis=1).flatten()
	    allnnout.append(nnout)
	    nnoutname.append("nnout_%s_M%d"%(key,mass))
    
    #for i in range(len(modellist)):
    for i in range(len(allnnout)):
	    #recfunctions.append_fields(base, names, data=None, dtypes=None, fill_value=-1, usemask=True, asrecarray=False)
	    a = recfunctions.append_fields(a, nnoutname[i], allnnout[i])
	    #a = recfunctions.(a, nnoutname[i], allnnout[i])

    #dataset.dtype.names = FullInputs + ['weight', 'mass'] + nnoutname
    #print "print a type ", type(a)," shape ",a.shape, " a.dtype.names ",a.dtype.names, " a[:10,:] ", a[:10]
    print "\n a dytpenames ", a.dtype.names 

    tfile = ROOT.TFile(outFile, "RECREATE")
    #ch_new = array2tree(dataset,  FullInputs + ['mass'] + nnoutname + ['finalweight'], "evtreeHME_nn")
    ch_new = array2tree(a, "evtreeHME_nn")
    
    ch_new.GetCurrentFile().Write() 
    ch_new.GetCurrentFile().Close()

def expectedSensitivity(hist_s, hist_b):
    xbins = hist_s.GetNbinsX() 
    hist_out = hist_s.Clone()
    for i in range (1, xbins+1):
        nsignal = hist_s.Integral(i, xbins)
        nbg = hist_b.Integral(i, xbins)
        if nbg <= .0000001:
            hist_out.SetBinContent(i, 5.0)
            hist_out.SetBinError(i, 5.0)
            continue
        hist_out.SetBinContent(i, nsignal/sqrt(nbg))
        error = sqrt(nsignal/nbg + (1.0/4.0)*nsignal*nsignal/(nbg*nbg))
        hist_out.SetBinError(i, error)
    return hist_out


def expectedSensitivity_Louvain(hist_s, hist_b):
    xbins = hist_s.GetNbinsX() 
    hist_out = hist_s.Clone()
    truenbins = (xbins)/3
    #controlregion1 = 25-3
    #controlregion2 = 50-6+1
    controlregion1 = truenbins
    controlregion2 = truenbins*2+1
    #print "expectedSensitivity_Louvain, controlregion1 ",controlregion1, " ", hist_s.GetBinCenter(controlregion1)," controlregion2 ",controlregion2," ", hist_s.GetBinCenter(controlregion2)

    for i in range (1, xbins+1):
        if i<= controlregion1 or i>= controlregion2:
            hist_out.SetBinContent(i,0)
            hist_out.SetBinError(i,0)
            continue
        nsignal = hist_s.Integral(i, controlregion2-1)
        nbg = hist_b.Integral(i, controlregion2-1)
        #print "i ",i ," bincenter ", hist_s.GetBinCenter(i)-controlregion1*0.04, " nbg ",nbg, " nsignal ", nsignal, " s over sqrt(B) ",nsignal/sqrt(nbg)
        if nbg <= .0000001:
            hist_out.SetBinContent(i, 5.0)
            hist_out.SetBinError(i, 5.0)
            continue
        hist_out.SetBinContent(i, nsignal/sqrt(nbg))
        error = sqrt(nsignal/nbg + (1.0/4.0)*nsignal*nsignal/(nbg*nbg))
        hist_out.SetBinError(i, error)
    return hist_out


def bestSensitivityNNWP(ch_signal, ch_bg, todraw, nbins, cuts_signal, cuts_bg, mass, plotname, outputdir):
    histbg = ROOT.TH1F("histbg","histbg",nbins, 0.0, 1.0)
    ch_bg.Draw(todraw+">> histbg",cuts_bg)
    histsignal = ROOT.TH1F("histsignal","histsignal",nbins, 0.0, 1.0)
    ch_signal.Draw(todraw+">> histsignal",cuts_signal)
    
    histsb = ROOT.TH1F("histsb","histsb",nbins, 0.0, 1.0)
    
    bestWP = 0.0
    bestS = 0.0
    for i in range(1, nbins+1):
	nbg = histbg.Integral(i, nbins)*1.0
	nsignal = histsignal.Integral(i, nbins)*1.0

        #print "i, ",i ," nbg ",nbg, " nsignal ",nsignal 
	if (nbg <= 0.0):
		continue
        histsb.SetBinContent(i, nsignal/sqrt(nbg+nsignal))
	if nsignal/sqrt(nbg+nsignal) > bestS:
		bestS = nsignal/sqrt(nbg+nsignal)
		bestWP = histbg.GetXaxis().GetBinCenter(i)
    if bestS == 0.0:
	    print "error in calculating S/sqrt(B) ",bestS," nncut ",bestWP
    else:
	    print "nnout ",todraw, " nncut ",bestWP," Sensitivity ",bestS 
    #makeplot = True
    if bestS>0.0:
	ctmp = ROOT.TCanvas("NNWP","NNWP",800, 600)
        ctmp.Clear()
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(.0)
        pad1.Draw()
        pad1.cd()
	#print "histbg ",histbg.Integral()," histsignal ",histsignal.Integral()
        maxy_s = histsignal.GetBinContent(histsignal.GetMaximumBin())
        maxy_b = histbg.GetBinContent(histsignal.GetMaximumBin())
        maxy = max(maxy_s, maxy_b)
	ratio = histsignal.Integral()/histbg.Integral()
	histbg.Scale(ratio)
        #histbg.Rebin(5)
        #histsignal.Rebin(5)
	histbg.SetLineColor(600)
	histbg.SetLineWidth(2)
	#histsignal.Scale(1.0/histsignal.Integral())
	histsignal.SetLineColor(632)
	histsignal.SetLineWidth(2)
	histbg.SetStats(0)
	histsignal.SetStats(0)
	histbg.Draw()
        histbg.SetMaximum(maxy*1.05)
	histsignal.Draw("same")
	histbg.SetTitle("Signal and Background DNN output distribution, Mass=%d GeV"%mass)
        histbg.GetYaxis().SetTitle("arbitary unit")
        #histbg.GetXaxis().SetTitle("DNN output")
	
        legend = ROOT.TLegend(0.35,0.63,0.70,0.88)
	legend.SetHeader("DNN training: kinematics+"+todraw[6:])
        legend.SetTextSize(0.04)
	ens = legend.AddEntry(histsignal, "Signal","pl")
	ens.SetTextColor(632)
	enbg = legend.AddEntry(histbg, "Background","pl")
	enbg.SetTextColor(600)
	legend.Draw("same")
        tex = ROOT.TLatex(0.35,0.55, "Maximum #frac{S}{#sqrt{B}} = %.2f @ NN = %.3f"%(bestS, bestWP))
        tex.SetNDC()
        tex.SetTextSize(.045)
        tex.Draw("same")

        ctmp.cd()
        ctmp.Update()
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, .29)
        pad2.SetTopMargin(0.)
        pad2.SetBottomMargin(.35)
        pad2.SetTicks(1,1)#draw x, y axis on both side (left right for y, and bottom up for x)
        pad2.Draw()
        pad2.cd()
        histsb.Draw("ep")
        histsb.SetStats(0)
        histsb.SetTitle("")

        histsb.GetXaxis().SetTitle("DNN output")
        histsb.GetXaxis().SetTitleSize(20)
        histsb.GetXaxis().SetTitleFont(43)
        histsb.GetXaxis().SetTitleOffset(3.0)
        histsb.GetXaxis().SetLabelSize(15)
        histsb.GetXaxis().SetLabelFont(43)#Absolute font size in pixel (precision 3)
        histsb.GetYaxis().SetTitle("S/#sqrt{B}")
        histsb.GetYaxis().SetNdivisions(505)
        histsb.GetYaxis().CenterTitle()
        histsb.GetYaxis().SetTitleSize(20)
        histsb.GetYaxis().SetTitleFont(43)
        histsb.GetYaxis().SetTitleOffset(.9)
        histsb.GetYaxis().SetLabelSize(15)
        histsb.GetYaxis().SetLabelFont(43)#Absolute font size in pixel (precision 3)

	#c1.SaveAs(outputdir+"bestSensitivityNNWP_%s_nnoutcut%d_bestS%d.png"%(todraw, bestWP*100, bestS*1000))
        #ctmp.SaveAs(outputdir+"bestSensitivityNNWP_%s_%s_nnoutcut%d.png"%(todraw, plotname, bestWP*100))
	ctmp.SaveAs(outputdir+"bestSensitivityNNWP_%s_%s_nnoutcut%d.pdf"%(todraw, plotname, bestWP*100))
    return bestWP, bestS

def hist1D(ch_signal, ch_bg, todraw, bins, xtitle, cuts_signallist, cuts_bglist, colors, title, legs, isNormalized, plotname):

    #color = [ROOT.kRed-4, ROOT.kBlue-4, ROOT.kMagenta-2, ROOT.kOrange-3, ROOT.kSpring, ROOT.kCyan, ROOT.kGreen+2]
    #color = [ROOT.kRed, ROOT.kBlue, ROOT.kMagenta+2, ROOT.kGreen+2, ROOT.kCyan+2, ROOT.kOrange+2, ROOT.kSpring]
    
    histbglist = []; histsignallist = []
    for i in range(len(cuts_signallist)):
	
        histbg = ROOT.TH1F("histbg%d"%i,"histbg",bins[0], bins[1], bins[2])
        ch_bg.Draw(todraw+">> "+histbg.GetName(),cuts_bglist[i])
        histsignal = ROOT.TH1F("histsignal%d"%i,"histsignal", bins[0], bins[1], bins[2])
        ch_signal.Draw(todraw+">> "+histsignal.GetName(),cuts_signallist[i])
	if histbg.GetEntries() == 0 or histsignal.GetEntries() == 0:
		print "error in ",plotname," todraw ",todraw," histbg entry ",histbg.GetEntries()," cut ",cuts_bglist[i], " hist_s entry ",histsignal.GetEntries(), " cut ",cuts_signallist[i] 
        if isNormalized:
	    histbg.Scale(1.0/histbg.Integral())
	    histsignal.Scale(1.0/histsignal.Integral())
	histbglist.append(histbg)
	histsignallist.append(histsignal)
    c1 = ROOT.TCanvas()
    c1.SetGridx(); c1.SetGridy(); c1.SetTickx(); c1.SetTicky()
    hs = ROOT.THStack("hs", title)

    deltay = 0.04*(len(cuts_signallist)+1)

    legend = ROOT.TLegend(0.45,0.65,0.88,0.66+deltay);
    legend.SetHeader("DNN training: kinematics+")
    #legend.SetFillColor(ROOT.kWhite); 
    #legend.SetTextSize(0.05); legend.SetTextFont(62)
    for i in range(len(cuts_signallist)):
	hs.Add(histbglist[i]) 
	hs.Add(histsignallist[i])
	histbglist[i].SetLineColor(colors[i])
	histbglist[i].SetLineStyle(1)
	histbglist[i].SetLineWidth(2)
	histsignallist[i].SetLineColor(colors[i])
	histsignallist[i].SetLineStyle(4)
	histsignallist[i].SetLineWidth(2)
	entry = legend.AddEntry(histsignallist[i], legs[i],"")
	entry.SetTextColor(colors[i])
	#legend.AddEntry(histbglist[i], legs[i]+",Background","pl")
    tex = ROOT.TLatex(0.6,0.2, "signal: - - - background: ---")
    tex.SetNDC()
    tex.SetTextSize(.035)
    hs.Draw("nostackhist")
    hs.GetHistogram().GetXaxis().SetTitle(xtitle)
    if isNormalized:
	    hs.GetHistogram().GetYaxis().SetTitle("Normalized to unity")

    legend.Draw("same")
    tex.Draw("same")
    c1.SaveAs(plotname+".C")
    c1.SaveAs(plotname+".pdf")

def fixhist2D(hist):
    xnbins = hist.GetNbinsX() 
    ynbins = hist.GetNbinsY() 
    for i in range(xnbins):
        for j in range(ynbins):
            if hist.GetBinContent(i, j)<0.0:
                hist.SetBinContent(i, j, 0.0)
def hist2D(ch, todraw, xbins, ybins, cuts, title, xtitle, ytitle, tex, plotname):
    c1 = ROOT.TCanvas()
    c1.Clear()
    hist = ROOT.TH2F("hist","hist",xbins[0], xbins[1], xbins[2], ybins[0], ybins[1], ybins[2])
    ch.Draw(todraw+">> hist", cuts)
    #print "todraw ",todraw, " cuts ",cuts," Integral ",hist.Integral()
    #warning: by using data driven method in DY estimation, some bin content after weighting is negative
    hist.SetMinimum(0.0)
    hist.Draw("colz")
    #hist.Draw()
    tex.Draw("same")
    hist.SetTitle(title)
    hist.GetXaxis().SetTitle(xtitle)
    hist.GetYaxis().SetTitle(ytitle)
    c1.SaveAs(plotname+".C")
    c1.SaveAs(plotname+".pdf")

def plotSB(file_slist, file_blist, nnlist, colors, title, mass, outdir, plotsuffix):
    ch_s = ROOT.TChain("evtreeHME_nn")
    ch_b = ROOT.TChain("evtreeHME_nn")
	
    for file_s in file_slist:
	ch_s.Add(file_s)
    for file_b in file_blist:
        ch_b.Add(file_b)
    print "file_s ",file_slist," file_b ",file_blist
    #bgweight = 708040.0/441719.0 ###  some jobs failed in processing ttbar 
    bgweight = 1.0#no job failed in this ntuples production
    variables = {"mt2":{
	    		"latexname":"M_{T2}",
			"bins":[80, 0.0, 500.0]
			},
		 "hme_h2mass_reco":{
			 "latexname":"HME Reco Mass",
			 "bins":[80, 200.0, 1000.0]
			 },
		 }
    legs = []
    nncutlist = []
    cuts_signallist = []
    cuts_backgroundlist = []
    #cuts_signal = "(35.87*1000*cross_section/event_weight_sum)*final_total_weight*(isMuMu && hme_h2mass_reco>250)"
    #cuts_bg = "(35.87*1000*cross_section/event_weight_sum)*final_total_weight*(isMuMu && hme_h2mass_reco>250)*%f"%(bgweight)
    signalweight = "final_total_weight"
    bgweight = "(35.87*1000*cross_section/event_weight_sum)*final_total_weight*%f"%(bgweight)
    for nnout in nnlist:
	cut = "isMuMu"# channel cut
	if "HME" in nnout:
            cut = cut+" && hme_h2mass_reco>250.0 "
	cuts_signal = "(%s)*(%s)"%(signalweight, cut)
	cuts_bg = "(%s)*(%s)"%(bgweight, cut)
	#cuts_signal = "(35.87*1000*cross_section/event_weight_sum)*final_total_weight*(isMuMu && hme_h2mass_reco>250)"
	#cuts_bg = "(35.87*1000*cross_section/event_weight_sum)*final_total_weight*(isMuMu && hme_h2mass_reco>250)*%f"%(bgweight)
	
	nncut, nnSensitivity = bestSensitivityNNWP(ch_s, ch_b, nnout, 200, cuts_signal, cuts_bg, mass, plotsuffix, outdir)
        legs.append(nnout[6:]+", DNN cut: %.3f"%(nncut))
	print "Singal ",plotsuffix," nncut ",nncut, " sensitivity ",nnSensitivity
	nncutlist.append(nncut)
	cuts_signallist.append("(%s)*(%s && %s>%f)"%(signalweight, cut, nnout, nncut))
	cuts_backgroundlist.append("(%s)*(%s && %s>%f)"%(bgweight, cut, nnout, nncut))
	#cuts_signallist.append("(35.87*1000*cross_section/event_weight_sum)*final_total_weight*(isMuMu && hme_h2mass_reco>250 && %s>%f)"%(nnout, nncut))
	#cuts_backgroundlist.append("(35.87*1000*cross_section/event_weight_sum)*final_total_weight*%f*(isMuMu && hme_h2mass_reco>250 && %s>%f)"%(bgweight, nnout, nncut))
    for var in variables:
	hist1D(ch_s, ch_b, var, variables[var]["bins"],  variables[var]["latexname"], cuts_signallist, cuts_backgroundlist, colors, title+","+variables[var]["latexname"], legs, True, os.path.join(outdir, var+"_afterDNN_"+plotsuffix))
	for i, nnout in enumerate(nnlist):
	    hist1D(ch_s, ch_b, var, variables[var]["bins"],  variables[var]["latexname"], [cuts_signallist[i]], [cuts_backgroundlist[i]], [colors[i]], title+","+variables[var]["latexname"], [legs[i]], True, os.path.join(outdir, var+"_afterDNN_"+plotsuffix+"_"+nnout))
    #### 2D hist
    for i, nnout in enumerate(nnlist):
	tex = ROOT.TLatex(0.7,0.2, nnout[6:])
        tex.SetNDC()
        tex.SetTextSize(.04)
        tex.SetTextFont(62)
	tex.SetTextColor(colors[i])
	hist2D(ch_s, "hme_h2mass_reco:mt2", variables["mt2"]["bins"], variables["hme_h2mass_reco"]["bins"], cuts_signallist[i]+"*(hme_h2mass_reco>0)", title+", plotting Signal", variables["mt2"]["latexname"], variables["hme_h2mass_reco"]["latexname"], tex,  os.path.join(outdir,"mt2VsHME_afterDNN_"+plotsuffix+"_"+nnout+"_signal"))
	hist2D(ch_b, "hme_h2mass_reco:mt2", variables["mt2"]["bins"], variables["hme_h2mass_reco"]["bins"], cuts_backgroundlist[i]+"*(hme_h2mass_reco>0)", title+", plotting Background", variables["mt2"]["latexname"], variables["hme_h2mass_reco"]["latexname"], tex, os.path.join(outdir,"mt2VsHME_afterDNN_"+plotsuffix+"_"+nnout+"_background"))
    for var in variables:
        for i, nnout in enumerate(nnlist):
	    tex = ROOT.TLatex(0.6,0.13, nnout[6:])
            tex.SetNDC()
            tex.SetTextSize(.04)
            tex.SetTextFont(62)
	    tex.SetTextColor(colors[i])

	    hist2D(ch_s, var+":"+nnout, [53, -0.02,1.04], variables[var]["bins"], signalweight+"*(%s>0.0)"%(var), title+", plotting Signal", "DNN output", variables[var]["latexname"], tex,  os.path.join(outdir, "nnoutVs"+var+"_beforeDNNcut_"+plotsuffix+"_"+nnout+"_signal"))
	    hist2D(ch_b, var+":"+nnout, [53, -0.02,1.04], variables[var]["bins"], bgweight+"*(%s>0.0)"%(var), title+", plotting Background", "DNN output", variables[var]["latexname"], tex,  os.path.join(outdir, "nnoutVs"+var+"_beforeDNNcut_"+plotsuffix+"_"+nnout+"_background"))


#extraweight is used for compensating the "HLT" safe ID cut, for Louvain ntuples
channelcuts = {"MuMu":{"cut":"isMuMu","extraweight": 24827.9/23891.14, "latex":"#mu#mu"},
		"MuEl":{"cut":"(isMuEl || isElMu)", "extraweight": 25786.8/30211.9, "latex":"e#mu"},
		"ElEl":{"cut":"isElEl","extraweight":8487.6/12190.1, "latex":"ee"},
		}
###
def histForlimits1D(backgrounddict, bgnames, filesignal, mass, todraw, cut, xbins, xtitle, suffix, outfile, plotname):

    if len(xbins) == 3:
        nbins = xbins[0]; xmin = xbins[1]; xmax =  xbins[2]
        xbins = []
        binwidth = (xmax-xmin)/nbins
        for i in range(0, nbins+1):
            xbins.append(xmin + i*binwidth)
        xbins = np.asarray(xbins)
    nnout = outfile.split("_")[-1][:-5]
    LouvainPlot = False
    #if (not LouvainPlot):
    #    nnout = ''
        
    if nnout == "MJJ":
        nnout = "MTandMT2MJJ"
    #print "nnout ",nnout
    chlist = {}
    #colors = [628, 596, 820, 432, 418]
    colors = [900-4, 820-3, 800-3, 860-3, 616-7, 432+2]
    treename = "evtreeHME_nn"
    #treename = "t"
    for key in backgrounddict:
	ch = ROOT.TChain(treename)
	for f in backgrounddict[key]:
	    ch.Add(f)
	chlist[key] = ch
    ch_s =  ROOT.TChain(treename); ch_s.Add(filesignal)
    rfile = ROOT.TFile(outfile, "UPDATE")
    xsec_signal = 5.0#pb
    weight = "(35.87*1000*cross_section/event_weight_sum)*final_total_weight"
    weight_s = "(35.87*1000*{xsec}/event_weight_sum)*final_total_weight".format(xsec = xsec_signal)
    for channel in channelcuts:
	allhists = []
        hs = ROOT.THStack("allbg_"+channel+"_"+suffix, "  ")
        legend = ROOT.TLegend(0.7,0.65,0.88,0.68+len(backgrounddict)*.06); 
	legend.SetTextSize(0.045); legend.SetTextFont(62)
        #legend.SetHeader("DNN training: kinematics+")
	BGSum = 0.0
        maxbgbin = 0.0
        #hist_data_fake = ROOT.TH1F("data_obs_"+channel+"_M%d_%s"%(mass, suffix), "data_obs_"+channel+"_M%d_%s"%(mass, suffix), len(xbins)-1, xbins)
	hist_data_fake = ROOT.TH1F("data_obs_"+channel+"_%s"%(suffix), "data_obs_"+channel+"_%s"%( suffix), len(xbins)-1, xbins)
	for i, key in enumerate(bgnames):
	    hist = ROOT.TH1F(key+"_"+channel+"_%s"%(suffix), key+"_"+channel+"_%s"%(suffix), len(xbins)-1, xbins)
	    finalcut = "("+ cut + " && "+ channelcuts[channel]["cut"] +")*"+weight+"*%f"%(channelcuts[channel]["extraweight"])
            print "todraw ",todraw," finalcut ",finalcut
	    chlist[key].Draw(todraw + ">> " + hist.GetName(), finalcut)
	    hs.Add(hist)
	    hist.SetFillColor(colors[i])
	    allhists.append(hist)
	    print "mass ",mass, " channel ", channel," bg ",key," rate ",hist.Integral()
            #if key == "TT":
            #    maxbgbin = hist.GetBinContent(hist.GetMaximumBin())
	    BGSum = BGSum + hist.Integral()
	    hist_data_fake.Add(hist)
	#for hist,name in zip(reversed(allhists), reversed(bgnames)):
	for i in range(len(allhists)):
	    legend.AddEntry(allhists[len(allhists)-i-1], bgnames[len(allhists)-i-1], "f")
	    allhists[len(allhists)-i-1].Write()
	hs.Write()
	
        hist_s = ROOT.TH1F("signal_"+channel+"_%s"%(suffix), "signal_"+channel+"_%s"%( suffix), len(xbins)-1, xbins)
	finalcut_s  = "("+ cut + " && "+ channelcuts[channel]["cut"] +")*"+weight_s+"*%f"%(channelcuts[channel]["extraweight"])
	ch_s.Draw(todraw + ">> " + hist_s.GetName(), finalcut_s)
        #hist_data_fake.Add(hist_s)
        maxbgbin = hist_data_fake.GetBinContent(hist_data_fake.GetMaximumBin())
        maxsignalbin = hist_s.GetBinContent(hist_s.GetMaximumBin())
        if maxsignalbin>maxbgbin:
            hs.SetMaximum(maxsignalbin*1.4)
        else:
            hs.SetMaximum(maxbgbin*1.3)

	print "mass ",mass, " channel ", channel," rate:signal ",hist_s.Integral()," BG(TT,DY,sT) ",BGSum," data(fake) ",hist_data_fake.Integral()

	hist_s.SetLineColor(1)
	hist_s.SetLineWidth(2)
	hist_s.Write()
	hist_data_fake.Write()
        if ( LouvainPlot):
            hist_sensitivity = expectedSensitivity_Louvain(hist_s, hist_data_fake)
        else:
            hist_sensitivity = expectedSensitivity(hist_s, hist_data_fake)

        hist_sensitivity.SetName("Sensitivity_"+channel+"_M%d"%mass)
        ens = legend.AddEntry(hist_s, "#splitline{Signal}{M=%d, %dpb}"%(mass, xsec_signal), "l")
	ens.SetTextSize(.03)

        c1 = ROOT.TCanvas("c", "canvas", 800, 800)
        c1.Clear()
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(.0)
        pad1.Draw()
        pad1.cd()
	hs.Draw("hist")
	hist_s.Draw("samehist")
	
        #tex1 = ROOT.TLatex(0.17,0.8, channelcuts[channel]["latex"]+" channel, "+nnout)
	tex1 = ROOT.TLatex(0.17,0.8, channelcuts[channel]["latex"]+" channel ")
	tex1.SetNDC(); tex1.SetTextSize(.045)
	tex2 = ROOT.TLatex(0.19,0.6, "M_{jj}<75 GeV "+"  "*12+" 75 GeV <=M_{jj}<140 GeV"+ "  "*12+" M_{jj} >= 140 GeV")
	tex2.SetNDC(); tex2.SetTextSize(.043)
	tex1.Draw("same")
        if ( LouvainPlot):
           tex2.Draw("same")
	#hs.GetHistogram().SetTitle(" #scale[1.4]{#font[61]{CMS}} Simulation Preliminary"+"  "*38+"35.87 fb^{-1} (13 TeV),2016")
        #hs.GetHistogram().SetTitle("")
	#hs.GetHistogram().SetTitleSize(.04)
	#hs.GetHistogram().SetTitleOffset(1.2)
	tex0 = ROOT.TLatex(0.1,0.92, " #scale[1.4]{#font[61]{CMS}} Simulation Preliminary"+"  "*24+"35.87 fb^{-1} (13 TeV),2016")
	tex0.SetNDC(); tex0.SetTextSize(.05); tex0.SetTextFont(42)
	tex0.Draw("same")
        #hs.GetHistogram().GetXaxis().SetTitle("DNN output, M_{jj} bins")
        #hs.GetHistogram().GetXaxis().SetTitle(xtitle)
        hs.GetHistogram().GetYaxis().SetTitle("Events")
        legend.Draw("same")
	#tex.Draw("same")

        c1.cd()
        c1.Update()
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, .29)
        pad2.SetTopMargin(0.)
        pad2.SetBottomMargin(.35)
        pad2.SetTicks(1,1)#draw x, y axis on both side (left right for y, and bottom up for x)
        pad2.Draw()
        pad2.cd()
        hist_sensitivity.Draw("ep")
        hist_sensitivity.SetStats(0)
        hist_sensitivity.SetTitle("")

        hist_sensitivity.GetXaxis().SetTitle(xtitle)
        hist_sensitivity.GetXaxis().SetTitleSize(20)
        hist_sensitivity.GetXaxis().SetTitleFont(43)
        hist_sensitivity.GetXaxis().SetTitleOffset(3.0)
        hist_sensitivity.GetXaxis().SetLabelSize(15)
        hist_sensitivity.GetXaxis().SetLabelFont(43)#Absolute font size in pixel (precision 3)
        hist_sensitivity.GetYaxis().SetTitle("S/#sqrt{B}")
        hist_sensitivity.GetYaxis().SetNdivisions(505)
        hist_sensitivity.GetYaxis().CenterTitle()
        hist_sensitivity.GetYaxis().SetTitleSize(20)
        hist_sensitivity.GetYaxis().SetTitleFont(43)
        hist_sensitivity.GetYaxis().SetTitleOffset(.9)
        hist_sensitivity.GetYaxis().SetLabelSize(15)
        hist_sensitivity.GetYaxis().SetLabelFont(43)#Absolute font size in pixel (precision 3)

        bestbin = hist_sensitivity.GetMaximumBin()
        bestWP = hist_sensitivity.GetBinCenter(bestbin)
        bestS = hist_sensitivity.GetBinContent(bestbin)
        tex_pad2 = ROOT.TLatex(0.2,0.35, "Maximum #frac{S}{#sqrt{B}} = %.1f @ %.3f"%(bestS, bestWP))
        tex_pad2.SetNDC()
        tex_pad2.SetTextSize(.035)
        #tex_pad2.Draw("same")

        c1.SaveAs(plotname+"_Radion_"+channel+"_"+suffix+".C")
        #c1.SaveAs(plotname+"_Radion_"+channel+".png")
        c1.SaveAs(plotname+"_Radion_"+channel+"_"+suffix+".pdf")
        #rfile.Close()
        print "done with histForlimits1D"
	

###
def histForlimits2D(backgrounddict, bgnames, filesignal, mass, xtodraw, ytodraw, cut, xbins, ybins, xtitle, ytitle, outfile, plotname):

    if len(xbins) == 3:
        nbins = xbins[0]; xmin = xbins[1]; xmax =  xbins[2]
        xbins = []
        binwidth = (xmax-xmin)*1.0/nbins
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

    chlist = {}
    #colors = [628, 596, 820, 432, 418]
    colors = [900-4, 820-3, 800-3, 860-3, 616-7, 432+2]
    treename = "evtreeHME_nn"
    for key in backgrounddict:
	ch = ROOT.TChain(treename)
	for f in backgrounddict[key]:
	    ch.Add(f)
	chlist[key] = ch
    ch_s =  ROOT.TChain(treename); ch_s.Add(filesignal)
    rfile = ROOT.TFile(outfile, "UPDATE")
    xsec_signal = 0.001
    weight = "(35.87*1000*cross_section/event_weight_sum)*final_total_weight"
    weight_s = "(35.87*1000*{xsec}/event_weight_sum)*final_total_weight".format(xsec = xsec_signal)
    todraw = ytodraw +":"+ xtodraw
    hist_0 = ROOT.TH2F("hist_0","", len(xbins)-1, xbins, len(ybins)-1, ybins)
    fixhist2D(hist_0)
    
    for channel in channelcuts:
	allhists = []
        hs = ROOT.THStack("allbg_M%d"%mass+channel, "  ")
        #legend = ROOT.TLegend(0.4,0.65,0.58,0.65+len(backgrounddict)*.06); 
        legend = ROOT.TLegend(0.4,0.7,0.5,0.75+len(backgrounddict)*.06); 
	legend.SetTextSize(0.045); legend.SetTextFont(62)
        #legend.SetHeader("DNN training: kinematics+")
	BGSum = 0.0
        maxbgbin = 0.0
	hist_data_fake = ROOT.TH2F("data_obs_"+channel+"_M%d"%mass, "data_obs_"+channel+"_M%d"%mass,  len(xbins)-1, xbins, len(ybins)-1, ybins)
	for i, key in enumerate(bgnames):
	    hist = ROOT.TH2F(key+"_"+channel+"_M%d"%mass, key+"_"+channel+"_M%d"%mass,  len(xbins)-1, xbins, len(ybins)-1, ybins)
	    finalcut = "("+ cut + " && "+ channelcuts[channel]["cut"] +")*"+weight+"*%f"%(channelcuts[channel]["extraweight"])
            print "todraw ",todraw," finalcut ",finalcut
	    chlist[key].Draw(todraw + ">> " + hist.GetName(), finalcut)
            if key == "DY":
                fixhist2D(hist)
            if key == "TT":
                maxbgbin = hist.GetBinContent(hist.GetMaximumBin())*1.05
	    hs.Add(hist)
	    hist.SetFillColor(colors[i])
	    allhists.append(hist)
	    print "mass ",mass, " channel ", channel," bg ",key," rate ",hist.Integral()
	    BGSum = BGSum + hist.Integral()
	    hist_data_fake.Add(hist)
	#for hist,name in zip(reversed(allhists), reversed(bgnames)):
	for i in range(len(allhists)):
	    legend.AddEntry(allhists[len(allhists)-i-1], bgnames[len(allhists)-i-1], "f")
	    allhists[len(allhists)-i-1].Write()
	hs.Write()
        
	
        hist_s = ROOT.TH2F("signal_"+channel+"_M%d"%mass, "signal_"+channel+"_M%d"%mass,  len(xbins)-1, xbins, len(ybins)-1, ybins)
	finalcut_s  = "("+ cut + " && "+ channelcuts[channel]["cut"] +")*"+weight_s+"*%f"%(channelcuts[channel]["extraweight"])
	ch_s.Draw(todraw + ">> " + hist_s.GetName(), finalcut_s)
        #hist_data_fake.Add(hist_s)
        maxsignalbin = hist_s.GetBinContent(hist_s.GetMaximumBin())

        if maxsignalbin>maxbgbin:
            hs.SetMaximum(maxsignalbin*1.4)
        else:
            hs.SetMaximum(maxbgbin*1.3)
        hs.SetMinimum(0.0)

	print "mass ",mass, " channel ", channel," rate:signal ",hist_s.Integral()," BG(TT,DY,sT) ",BGSum," data(fake) ",hist_data_fake.Integral()
        if hist_data_fake.GetEntries() <= 0:
            print "hist_data_fake entries is not greater than 0!!! error !!"
            raise ValueError
        #hist_s.SetLineColor(400)
        hist_s.SetFillColor(612)
	hist_s.SetLineWidth(2)
        hist_s.SetStats(0)
        hist_s.GetXaxis().SetTitle(xtitle)
        hist_s.GetXaxis().SetTitleOffset(1.5)
        hist_s.GetXaxis().CenterTitle()
        hist_s.GetYaxis().SetTitle(ytitle)
        hist_s.GetYaxis().SetTitleOffset(1.5)
        hist_s.GetYaxis().CenterTitle()
        hist_data_fake.SetStats(0)
        hist_data_fake.GetXaxis().SetTitle(xtitle)
        hist_data_fake.GetXaxis().SetTitleOffset(1.5)
        hist_data_fake.GetXaxis().CenterTitle()
        hist_data_fake.GetYaxis().SetTitle(ytitle)
        hist_data_fake.GetYaxis().SetTitleOffset(1.5)
        hist_data_fake.GetYaxis().CenterTitle()
	hist_s.Write()
	hist_data_fake.Write()
        ens = legend.AddEntry(hist_s, "#splitline{Signal}{M=%d, %dpb}"%(mass, xsec_signal), "f")
	ens.SetTextSize(.03)

        c1 = ROOT.TCanvas("c1","c1",800, 400)
        c1.Clear()
        c1.Divide(2)
        c1.cd(1)
	hs.Draw("LEGO1 0")
	
	tex1 = ROOT.TLatex(0.17,0.85, channelcuts[channel]["latex"]+" channel")
	tex1.SetNDC(); tex1.SetTextSize(.045)
	tex2 = ROOT.TLatex(0.19,0.5, "M_{jj}<75 GeV "+"  "*5+" 75 GeV <=M_{jj}<140 GeV"+ "  "*5+" M_{jj} >= 140 GeV")
	tex2.SetNDC(); tex2.SetTextSize(.033)
        #tex1.Draw("same")
        #tex2.Draw("same")
	#hs.GetHistogram().SetTitle(" #scale[1.4]{#font[61]{CMS}} Simulation Preliminary"+"  "*38+"35.87 fb^{-1} (13 TeV),2016")
        channelname = channelcuts[channel]["latex"]+" channel"
	hs.SetTitle("Background, "+channelname)
	#hs.GetHistogram().SetTitleSize(.04)
	#hs.GetHistogram().SetTitleOffset(1.2)
	tex0 = ROOT.TLatex(0.1,0.91, " #scale[1.4]{#font[61]{CMS}} Simulation Preliminary"+"  "*14+"35.87 fb^{-1} (13 TeV),2016")
	tex0.SetNDC(); tex0.SetTextSize(.04); tex0.SetTextFont(42)
        #tex0.Draw("same")
        #hs.GetHistogram().GetXaxis().SetTitle("DNN output, M_{jj} bins")
        hs.GetHistogram().GetXaxis().SetTitle(xtitle)
        hs.GetHistogram().GetXaxis().SetTitleOffset(1.5)
        hs.GetHistogram().GetXaxis().CenterTitle()

	hs.GetHistogram().GetYaxis().SetTitle(ytitle)
        hs.GetHistogram().GetYaxis().SetTitleOffset(1.5)
        hs.GetHistogram().GetYaxis().CenterTitle()

        c1.cd(2)
        hs_s = ROOT.THStack("hs_s","Signal, %d GeV, "%(mass)+channelname)
        #hist_0.SetFillColor(596)
        hs_s.Add(hist_s)
        #hs_s.Add(hist_0)
        hs_s.Draw("LEGO1 0")
        hs_s.GetHistogram().GetXaxis().SetTitle(xtitle)
        hs_s.GetHistogram().GetXaxis().SetTitleOffset(1.5)
        hs_s.GetHistogram().GetXaxis().CenterTitle()

	hs_s.GetHistogram().GetYaxis().SetTitle(ytitle)
        hs_s.GetHistogram().GetYaxis().SetTitleOffset(1.5)
        hs_s.GetHistogram().GetYaxis().CenterTitle()
        c1.Update()
        c1.cd()

        legend.Draw("same")
	#tex.Draw("same")
        c1.SaveAs(plotname+"_RadionM%d_"%mass+channel+".C")
        c1.SaveAs(plotname+"_RadionM%d_"%mass+channel+".pdf")

        c2 = ROOT.TCanvas("c2","c2",800, 400)
        c2.Clear()
        c2.Divide(2,1)
        c2.cd(1)
        #hs.SetTitle("Background, "+channelname)
        hist_data_fake.SetMinimum(0.0)
        #hist_data_fake.SetMaximum()
        hist_data_fake.Draw("colz")
        hist_data_fake.SetTitle("Background, "+channelname)
        c2.cd(2)
        hist_s.Draw("colz")
        hist_s.SetTitle("Signal, %s GeV, "%(mass)+channelname)
        c2.SaveAs(plotname+"_RadionM%d_"%mass+channel+"_2plots.C")
        c2.SaveAs(plotname+"_RadionM%d_"%mass+channel+"_2plots.pdf")
        

	


#filepath = "/Users/taohuang/Documents/DiHiggs/20170530/20171021_Louvain/"
#filepath = "/Users/taohuang/Documents/DiHiggs/20170530/testfile/"
#filepath = "/Users/taohuang/Documents/DiHiggs/20170530/20171102_Louvain/"
#filepath = "/Users/taohuang/Documents/DiHiggs/20180126/20180205_20180202_10k_Louvain_ALL/"
filepath = "/Users/taohuang/Documents/DiHiggs/20180316_NanoAOD/HHNtuple_20180518_addSys/"
filepath = "/Users/taohuang/Documents/DiHiggs/20180316_NanoAOD/20180619_HHbbWW_addHME_addSys_10k/"
#outdir = "/Users/taohuang/Documents/DiHiggs/20170530/20171021_Louvain_addNN/"
outdir = filepath
suffix = 'addNN'
date_suffix = '{:%Y-%m-%d}'.format(datetime.date.today())
output_suffix = '2018-06-19_test_addNN'
output_folder = os.path.join(outdir, output_suffix)
output_folder = "/Users/taohuang/Documents/DiHiggs/20180316_NanoAOD/20180619_HHbbWW_addHME_addSys_10k_addNN/"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

modellist = ['MTonly','MTandMT2','MTandMT2_MJJ', 'MTandMT2_HME']
#modellist = ['MTonly','MT2only','MTandMT2','MTandMT2_MJJ','MTandMT2_HME','MTandMT2_HMEMJJ', 'MTandMJJ']
modellist = ['MTonly']
allfiles = os.listdir(filepath)
#allfiles = ['TT_all.root']
def makeNtuple_prediction(masses):
     
    for f in allfiles:
	if f.endswith(".root"):
            #name = f.split('.')[0][:-4]
            name = f.split('.')[0]
            print "file ",f," samplename ",name
	    outname = ''
	    try:
	    	sampletype =  SampleTypeLUT.sampletypedict[name]
	    except KeyError:
	    	##background
                #sampletype =  SampleTypeLUT.sampletypedict[name[:-4]]
	    	#outname = name+"_addNN_M%s.root"%mass
                #print "Failed to find sample type ",name
                sampletype = 1
            
	    #addsampleType(filepath+f, sampletype, name+"_addType.root")
	    #weight = "event_weight * trigeff * jjbtag_heavy * jjbtag_light * llidiso * pu"
            #weight = "total_weight"
            weight = "sample_weight*event_reco_weight"
            #if name.startswith("DY"):
            #	weight = weight+"* dy_nobtag_to_btagM_weight"
	    print 'file ',file,' name ',name,' type ',type(sampletype)," sampletype ",sampletype," weight ",weight
	    if sampletype>=100:
	        outname = name+"_addNN.root"
	        output_file = os.path.join(output_folder, outname)		
	        mass = SampleTypeLUT.getSignalMass(sampletype)
	        writeNNToTree(filepath+f, weight, modellist, SampleTypeLUT.ModelLUT, sampletype, [mass], output_file)
	    else:
	    	#for mass in masses:
	        #outname = name+"_addNN_M%s.root"%mass
	        outname = name+"_addNN.root"
	        output_file = os.path.join(output_folder, outname)		
	        writeNNToTree(filepath+f, weight, modellist, SampleTypeLUT.ModelLUT, sampletype, masses, output_file)


#makeNtuple_prediction([260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900])
#makeNtuple_prediction([260, 400])

#nnlist = ["nnout_MTonly", "nnout_MT2only","nnout_MTandMT2","nnout_MTandMT2_MJJ","nnout_MTandMT2_HME","nnout_MTandMT2_HMEMJJ","nnout_MTandMJJ"]
#nnlist = ["nnout_MTonly", "nnout_MTandMJJ"]
nnlist = ["nnout_MTonly"]
#plotdir = "plotsafterDNN/"
def makeSBplots(nnlist, plotdir, masspoints):
    os.system("mkdir -p "+plotdir)
    #masspoints = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 800, 900]
    #colors = [632, 600, 618, 418, 434, 802]
    #colors = [628, 596, 614, 797, 820, 432, 418]
    colors = [900-4, 820-3, 800-3, 860-3, 616-7, 432+2]
    for mass in masspoints: 
	file_slist = []
        file_slist.append(output_folder+"/radion_M%d_addNN.root"%(mass))
	bgnamelist = ["TT","DYM10to50","DYToLL0J","DYToLL1J","DYToLL2J","sT_top","sT_antitop"]
	#bgnamelist = ["DYM10to50","DYToLL0J","DYToLL1J","DYToLL2J"]
	file_blist = []
	for bgname in bgnamelist:
		file_blist.append( output_folder+"/%s_all_addNN.root"%(bgname) )
	title = "Signal, Radion M=%d GeV "%mass
	thisnnlist = []
	for nn in nnlist:
		thisnnlist.append(nn+"_M%d"%mass)
        plotSB(file_slist, file_blist, thisnnlist, colors, title, mass, plotdir, "radion_M%d"%mass)

#makeSBplots(nnlist, "plotsafterDNN_splitted_20200331_test/", [270, 750])
#makeSBplots(nnlist, "plotsafterDNN_splitted_20171102/", [260, 350, 400, 500])


def makeBackgroundshist(filefolder, masspoints, variable, nbins, xtitle, outdir):

    bgnames = ["sT","DY","TT"]
    #bgnames = ["TT"]
    outfile = os.path.join(outdir, "Backgrounds_signal_allinputs.root")
    plotname = os.path.join(outdir, "Kinematics_%s"%variable)
    ###create tfile
    tfile = ROOT.TFile(outfile, "RECREATE")
    tfile.Close()
    todraw = variable
    for mass in masspoints:
        cut = "1"
        #todraw = "(({nnout}_M{mass}-3.0/25)*(jj_M<75 && {nnout}_M{mass}>3.0/25)+(jj_M>=75 && jj_M<140 && {nnout}_M{mass}>3.0/25)*({nnout}_M{mass}+1-6.0/25)+(jj_M>=140 && {nnout}_M{mass}>3.0/25)*({nnout}_M{mass}+2-9.0/25))".format(nnout = nnout, mass=mass)
        file_s = os.path.join(filefolder, "radion_M%d_addNN.root"%(mass))
        file_s = os.path.join(filefolder, "GluGluToRadionToHHTo2B2VTo2L2Nu_M-%d_narrow_13TeV-madgraph-v2_HME_Friends_NN.root"%mass)
        backgrounddict = {}
        """
        for bgname in bgnames:
            backgrounddict[bgname] = []
	file_s = os.path.join(filefolder, "radion_M%d_addNN.root"%(mass))
        for background in backgrounds:
            bginput = background['input']
            files = get_files_from_glob(os.path.join(INPUT_FOLDER, background['input']))
            bgname = bginput[:2]
            print "bginput ",bginput," bgname ",bgname, " files ", files
            backgrounddict[bgname]  = backgrounddict[bgname]  + files
        for bgname in bgnames:
            print "bgname ",bgname, " files ",backgrounddict[bgname]

        """
	for bg in bgnames:
	    allfiles = os.listdir(filefolder)
	    backgrounddict[bg] = []
	    for f in allfiles:
		if f.startswith(bg) and f.endswith(".root"):
		   backgrounddict[bg].append(os.path.join(filefolder, f))
        suffix = ''
        histForlimits1D(backgrounddict, bgnames, file_s, mass, todraw, cut, nbins, xtitle, suffix, outfile, plotname)

variablesdir = "Background_signal_xsec5pb_MTonlyLouvainTraining_Test"
os.system("mkdir -p "+variablesdir)
varibales = ['jj_pt', 'll_pt', 'll_M', 'll_DR_l_l', 'jj_DR_j_j', 'llmetjj_DPhi_ll_jj', 'llmetjj_minDR_l_j', 'llmetjj_MTformula','mt2', 'jj_M','hme_h2mass_reco']
#variables = ['lep1_pt']
#makeBackgroundshist(output_folder, [400], 'llmetjj_MTformula', [50, 0.0, 500.0],"MT", variablesdir)
def plotallkinematics():
    #output_folder = "/Users/taohuang/Documents/DiHiggs/20180316_NanoAOD/HHNtuple_20180328_fixedleptonDZeff"
    print "Ntuple folder ",output_folder
    
    makeBackgroundshist(output_folder, [400], 'll_M', [50, 12.0, 76.0], "M_{ll}", variablesdir)
    #makeBackgroundshist(output_folder, [400], 'll_DR_l_l', [50, .0, 6.0], "#DeltaR_{ll}", variablesdir)
    #makeBackgroundshist(output_folder, [400], 'jj_M', [50, 0.0, 400.0], "M_{jj}",variablesdir)
    #makeBackgroundshist(output_folder, [400], 'jj_DR_j_j', [50, .0, 6.0], "#DeltaR_{jj}",variablesdir)
    #makeBackgroundshist(output_folder, [400], 'llmetjj_DPhi_ll_jj', [24, .0, 3.1415926],"#Delta#phi(ll,jj)", variablesdir)
    #makeBackgroundshist(output_folder, [400], 'll_pt', [50, 0.0, 450.0], "Dilepton p_{T}", variablesdir)
    #makeBackgroundshist(output_folder, [400], 'jj_pt', [50, 0.0, 450.0], "Dijet p_{T}", variablesdir)
    #makeBackgroundshist(output_folder, [400], 'llmetjj_minDR_l_j', [50, .0, 5.0], "#DeltaR_{l,j}", variablesdir)
    #makeBackgroundshist(output_folder, [400], 'llmetjj_MTformula', [50, 0.0, 500.0],"MT", variablesdir)
    #makeBackgroundshist(output_folder, [400], 'met_pt', [50, 0.0, 500.0],"MET", variablesdir)
#plotallkinematics()

def generateHMEbins(mass):
    lowM = 250.0; highM = 1200.0
    xbins = [lowM]
    x = lowM
    step = 25.0 
    if mass >= 400 and mass < 500:
        step = 30.0
    elif mass>=500 and mass < 700:
        step = 40.0
    elif mass>=700:
        step = 50.0
    gap1 = 50 + mass*.1 #120.0
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
    xbins.append(1200.0)
    print "Benchmark ", mass ," HME mass bins ",xbins
    return np.asarray(xbins)
            

def makeLouvainLimithist(masspoints, nnout, outdir):
    bgnames = ["sT","DY","TT"]
    outfile = os.path.join(outdir, "Hhh_FinalBGYield_xsec1pb_NN_%s_%s.root"%(date_suffix, nnout))
    plotname = os.path.join(outdir, "Hhh_FinalBGYield_xsec1pb_NN_%s"%nnout)
    ###create tfile
    tfile = ROOT.TFile(outfile, "RECREATE")
    tfile.Close()
    nnbins = []
    for x in range(3, 75-6+1):
        nnbins.append(x*1.0/25.0)
    print "nnbins ",nnbins
    nnbins_x = np.asarray(nnbins)
    xtitle = "DNN output, Mjj bin"
    for mass in masspoints:
        cut = "({nnout}_M{mass}>3.0/25 && hme_h2mass_reco>=250)".format(nnout = nnout, mass = mass)
        todraw = "(({nnout}_M{mass})*(jj_M<75 && {nnout}_M{mass}>3.0/25)+(jj_M>=75 && jj_M<140 && {nnout}_M{mass}>3.0/25)*({nnout}_M{mass}+1-3.0/25)+(jj_M>=140 && {nnout}_M{mass}>3.0/25)*({nnout}_M{mass}+2-6.0/25))".format(nnout = nnout, mass=mass)
	file_s = os.path.join(output_folder, "radion_M%d_addNN.root"%(mass))
        print "file_s ",file_s
        backgrounddict = {}
	for bg in bgnames:
	    allfiles = os.listdir(output_folder)
	    backgrounddict[bg] = []
	    for f in allfiles:
		if f.startswith(bg):
		   backgrounddict[bg].append(os.path.join(output_folder, f))
        suffix = "M%d"%(mass)
        histForlimits1D(backgrounddict, bgnames, file_s, mass, todraw, cut, nnbins_x, xtitle,suffix, outfile, plotname)
        


def make1DLimithist(masspoints, nnout, nnoutcut, xtodraw, suffix, outfile, outdir):
    bgnames = ["sT","DY","TT"]
    plotname = os.path.join(outdir, "Hhh_FinalBGYield_xsec1pb_HME_%s"%nnout)
    ###create tfile
    xtitle = "HME, Mjj bin"
    for mass in masspoints:
        #cut = "({nnout}_M{mass}>{nnoutcut}) && {xvariable}>250.0".format(nnout = nnout, mass = mass, nnoutcut = nnoutcut, xvariable = xtodraw)
        cut = "{xvariable}>250.0".format(xvariable = xtodraw)
        nncut = "&& ("
        for channel in nnoutcut:
            nncut = nncut+" ( {ch} && {nnout}_M{mass}>{nnoutcut}) || ".format(ch = channel, nnout = nnout, mass = mass, nnoutcut = nnoutcut[channel])
        cut = cut+ nncut[:-3] + ") "
        print "cut ", cut
        xbins = generateHMEbins(mass)
        xmin = xbins[0]; xmax = xbins[-1]; nxbin = len(xbins)
        todraw = "(({xvariable})*(jj_M<75)+(jj_M>=75 && jj_M<140)*({xvariable}+{massgap})+(jj_M>=140)*({xvariable}+2*{massgap}))".format(xvariable = xtodraw, massgap = xmax-xmin)
        tmpbins = xbins[1:]
        for i in range (1, 3):
            for xval in tmpbins:
                #print "xval + i*(xmax-xmin) ",xval + i*(xmax-xmin)," xval ",xval," xmax -xmin ",xmax-xmin
                xbins = np.append(xbins, xval + i*(xmax-xmin))

        #todraw = "((jj_M>=75 && jj_M<140)*({xvariable}))".format(xvariable = xtodraw)
        #print "final xbins ",xbins
        #file_s = os.path.join(output_folder, "radion_M%d_addNN.root"%(mass))
        file_s = os.path.join(output_folder, "GluGluToRadionToHHTo2B2VTo2L2Nu_M-%d_narrow_13TeV-madgraph-v2_HME_Friends_NN.root"%mass)
        backgrounddict = {}
	for bg in bgnames:
	    allfiles = os.listdir(output_folder)
	    backgrounddict[bg] = []
	    for f in allfiles:
		if f.startswith(bg):
		   backgrounddict[bg].append(os.path.join(output_folder, f))
        histForlimits1D(backgrounddict, bgnames, file_s, mass, todraw, cut, xbins, xtitle, suffix, outfile, plotname)

def make1DLimithist_HEMwithNNcut(masspoints, nnout, outdir):
    bgnames = ["sT","DY","TT"]
    treename = "evtreeHME_nn"
    bestNNcut = {}
    bestNNSoverB = {}

    outfile = os.path.join(outdir, "Hhh_FinalBGYield_xsec1pb_HMENNcut_%s_%s.root"%(date_suffix, nnout))
    tfile = ROOT.TFile(outfile, "RECREATE")
    tfile.Close()
    for mass in masspoints:
      bestNNcut[mass] = {}
      bestNNSoverB[mass] =  {}
      ch_b = ROOT.TChain(treename)
      backgrounddict = {}
      for bg in bgnames:
          allfiles = os.listdir(output_folder)
          for f in allfiles:
              if f.startswith(bg):
                  bgfile = os.path.join(output_folder, f)
                  ch_b.Add(bgfile)
      for channel in ["isMuMu", "isElMu","isElEl"]:
        ch_s = ROOT.TChain(treename)
        signalweight = "final_total_weight"
        jobsdone = 1.0#
        ### for each channel ????
        bgweight = "(35.87*1000*cross_section/event_weight_sum)*final_total_weight*%f"%(jobsdone)
        cut = "(jj_M>=75 && jj_M<140 && hme_h2mass_reco>=250 && %s )"%channel
        #cut = "(jj_M>=75 && jj_M<140) && %s"%channel
        xbins = "(25, 0, 1.0)"
        file_s =  os.path.join(output_folder, "radion_M%d_addNN.root"%(mass))
        ch_s.Add(file_s)
        cuts_signal = "(%s)*(%s)"%(cut, signalweight) 
        cuts_bg = "(%s)*(%s)"%(cut, bgweight) 
        #bestWP, bestS = bestSensitivityNNWP(ch_s, ch_b, nnout+"_M%d"%mass, 80, cuts_signal, cuts_bg, mass, channel, outdir)
        #bestNNcut[mass][channel] = bestWP
        #bestNNSoverB[mass][channel] = bestS

    #print "bestNNcut ",bestNNcut, " bestNNSoverB ",bestNNSoverB
    for mass in masspoints:
        outdir_tmp = os.path.join(outdir, "M_%d"%mass)
        if not os.path.exists(outdir_tmp):
            os.makedirs(outdir_tmp)
        bestWPs = bestNNcut[mass] 
        for WP in range(0, 7):
            step = .1*WP
            #nnoutcut = {"isMuMu":0.12+step, "isElMu":0.12+step, "isElEl":0.12+step}
            nnoutcut = {"isMuMu":step, "isElMu":step, "isElEl":step}
            suffix = "RadionM%d_WP%d"%(mass, WP)
            make1DLimithist([mass], nnout, nnoutcut, "hme_h2mass_reco", suffix,outfile, outdir_tmp)
        """
        nnoutcut = {}
        NWPs = 1
        for i in range(0, NWPs):
            for channel in bestWPs:
                nnoutcut[channel] = bestWPs[channel] - i*0.04
                if i == NWPs-1:
                    nnoutcut[channel] = 0.12
                if nnoutcut[channel] <= 0.0: 
                    print "warning, nnoutcut less than zero !!!, pass"
                    continue
            suffix = "RadionM%d_WP%d"%(mass, i)
            #make1DLimithist([mass], nnout, nnoutcut, "hme_h2mass_reco", suffix,outfile, outdir_tmp)


        """
def make2DLimithist(masspoints, nnout, ytodraw, ytitle, outdir):
    bgnames = ["sT","DY","TT"]
    outfile = os.path.join(outdir, "Hhh_FinalBGYield_2Dlimits_%s_%s.root"%(date_suffix, nnout))
    plotname = os.path.join(outdir, "Hhh_FinalBGYield_2Dlimits_%s_%s"%(nnout, ytodraw))
    ###create tfile
    tfile = ROOT.TFile(outfile, "RECREATE")
    tfile.Close()
    #nnbins = [75, 0.0, 3.0] 
    nnstep = 0.05
    nnbins = [int(3.0/nnstep),0.,3.0]
    nncut = 0.35
    nnbins[0] = nnbins[0] - 3*int(nncut/nnstep)
    nnbins[1] = nnbins[1]+nncut
    nnbins[2] = nnbins[2] - nncut*2
    print "nnbins ",nnbins, " nncut ",nncut
    for mass in masspoints:
        massbins = generateHMEbins(mass)
        cut = "({nnout}_M{mass} > {nncut}) && {y}>0".format(nnout = nnout, mass = mass, y=ytodraw, nncut = nncut)
        xtodraw  = "(({nnout}_M{mass})*(jj_M<75) + (jj_M>=75 && jj_M<140)*({nnout}_M{mass}+1-1*{nncut})+(jj_M>=140)*({nnout}_M{mass}+2- 2*{nncut}))".format(nnout = nnout, mass=mass, nncut = nncut)
        xtitle = "DNN output, Mjj bin"
        #file_s = os.path.join(output_folder, "radion_M%d_addNN.root"%(mass))
        #file_s = os.path.join(output_folder, "radion_M%d_addNN.root"%(mass))
        file_s = os.path.join(output_folder, "GluGluToRadionToHHTo2B2VTo2L2Nu_M-%d_narrow_13TeV-madgraph-v2_HME_Friends_NN.root"%mass)
        backgrounddict = {}
	for bg in bgnames:
	    allfiles = os.listdir(output_folder)
	    backgrounddict[bg] = []
	    for f in allfiles:
		if f.startswith(bg):
		   backgrounddict[bg].append(os.path.join(output_folder, f))
        print("backgrounddict  ", backgrounddict)
        histForlimits2D(backgrounddict, bgnames, file_s, mass, xtodraw,ytodraw, cut, nnbins, massbins, xtitle, ytitle, outfile, plotname)

def writeintoworkspace(masslist, infile, prefix_out, outdir):
    for mass in masslist:
        outfile = os.path.join(outdir, prefix_out+"_M%d_workspace.root"%mass)
        command = 'root -b -q "rf_workspace.C(%d, \\"%s\\", \\"%s\\")"'%(mass, infile, outfile)
        print "writeintoworkspace command ",command
        os.system(command)
#limitdir = "HistForLimits_2D_xsec1pb/"
#limitdir = "HistForLimits_HME_1D_xsec1pb_nnout5_MTonly/"
limitdir = "HistForLimits_1D_xsec5pb_DNNshape_2020Test/"
#limitdir = "MTandMT2_HME_2D_xsec1fb_20180312_2D_NNbin10_nncut035_test/"
#limitdir = "HistForLimits_1D_xsec1pb_HME_WP1_newbins/"
#os.system("mkdir -p "+limitdir)
#for mass in [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900]:
    #generateHMEbins(mass)
#generateHMEbins(800)
makeLouvainLimithist([260], "nnout_MTonly",limitdir)
#make1DLimithist([300], "nnout_MTonly", 0.2, "hme_h2mass_reco", limitdir)
#make2DLimithist([400], "nnout_MTonly", "hme_h2mass_reco", "HME reco mass", limitdir)
#make2DLimithist([260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900], "nnout_MTonly", "hme_h2mass_reco", [80, 200.0, 1000.0], "HME reco mass", limitdir)
model = "MTonly"
#plotmodel(SampleTypeLUT.ModelLUT, model, model+"_visualization.pdf")
#plotmodel(SampleTypeLUT.ModelLUT, model, model+"_visualization.png")

#for nnout in modellist:
#for nnout in ["MTandMT2","MTandMT2_MJJ"]:
for nnout in ["MTonly", "MTandMT2"]:
#for nnout in ["MTonly"]:
    #masspoints = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900]
    masspoints = [400]
    #masspoints = [260, 300, 350, 400, 500, 600, 750, 900]

    #makeLouvainLimithist(masspoints, "nnout_"+nnout, limitdir)

    #nnoutcut = 0.0
    #xtodraw = "hme_h2mass_reco"
    #make1DLimithist(masspoints, "nnout_"+nnout, nnoutcut, xtodraw, limitdir)
    #make1DLimithist_HEMwithNNcut(masspoints, "nnout_"+nnout, limitdir)
    #make2DLimithist(masspoints,  "nnout_"+nnout, "hme_h2mass_reco", "HME reco mass", limitdir)
    
    #infile = limitdir + "Hhh_FinalBGYield_2Dlimits_2018-03-13_nnout_%s.root"%(nnout)
    #writeintoworkspace(masspoints, infile, "HME_%s_2D_bin10"%(nnout), limitdir)

    print "nnout ",nnout

print "done with writeNNToTee.py"
