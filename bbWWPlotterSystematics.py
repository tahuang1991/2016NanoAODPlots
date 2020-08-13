import os
import sys
import ROOT
from math import fabs, sqrt

## ------------------------------ ##
## Setup Matplotlib Configuration ##



class bbWWPlotterSystematics(object):
    """Class for handling systematic uncertainties to put in plots"""
    def __init__(self,  inputfiledict, nominalWeight, systematiclist, DYdatadriven):
        """Initialize class for each sample you want systematics"""
        self.channels = ["MuMu", "MuEl","ElEl"]
        self.channelcuts = {"MuMu":"isMuMu", "MuEl":"(isMuEl || isElMu)","ElEl":"isElEl"}
        #self.outfile = outfile
        #self.addStatistics = True
        #self.writeSystematichists = False
        
        
        self.SysDict = {
                "CMS_eff_b_heavy":{  "weight":"jjbtag_heavy_nominal",
                                     "description":"b-tagging uncertainty from b,c-jets", 
                                     "event_weight_sum": "nominal",
                                     },
                "CMS_eff_b_light":{  "weight":"jjbtag_light_nominal",
                                     "description":"b-tagging uncertainty from light-jets",
                                     "event_weight_sum": "nominal",
                                     },
                "CMS_pu":{ "weight":"event_pu_weight_nominal",
                           "description":"pielup uncertainty",
                           "event_weight_sum":"nominal",
                           },
                "CMS_pdf":{"weight":"event_pdf_weight_nominal",
                            "description":"PDF uncertainty",
                            "event_weight_sum":"nominal",
                            },
                "CMS_eff_trigger":{"weight":"lep1trgsf_nominal * lep2trgsf_nominal", 
                                   "description": "trigger eff uncertainty",
                                   "event_weight_sum":"nominal",
                                   },
                #"CMS_eff_e":{"weight": "(isMuEl*lep2trackingsf_nominal*lep2HLTsafeIDsf_nominal*lep2IDsf_nominal + isElMu*lep1trackingsf_nominal*lep1HLTsafeIDsf_nominal*lep1IDsf_nominal+ isElEl*lep1trackingsf_nominal*lep1HLTsafeIDsf_nominal*lep1IDsf_nominal*lep2trackingsf_nominal*lep2HLTsafeIDsf_nominal*lep2IDsf_nominal+isMuMu)", 
                "CMS_eff_e":{"weight": "(isMuEl*lep2trackingsf_nominal*lep2IDsf_nominal + isElMu*lep1trackingsf_nominal*lep1IDsf_nominal+ isElEl*lep1trackingsf_nominal*lep1IDsf_nominal*lep2trackingsf_nominal*lep2IDsf_nominal+isMuMu)", 
                            "description":"Electron eff",
                            "event_weight_sum":"nominal",
                            },
                "CMS_eff_mu":{"weight": "(isElMu*lep2trackingsf_nominal*lep2IDsf_nominal + isMuEl*lep1trackingsf_nominal*lep1IDsf_nominal+ isMuMu*lep1trackingsf_nominal*lep1IDsf_nominal*lep2trackingsf_nominal*lep2IDsf_nominal+isElEl)",
                            "description":"Muon eff",
                            },
                "CMS_iso_mu":{"weight": "(lep1Isosf_nominal * lep2Isosf_nominal)", 
                              "description": "Muon Isolation",
                              "event_weight_sum":"nominal",
                              },
                "QCDscale":{"weight":None,
                            "description":"QCD scale uncertainty",
                            "event_weight_sum":None,
                            },
                "MC_statistical":{
                    "weight": None,
                    "description":"uncertainty due to finite MC statistics",
                    "event_weight_sum":"nominal",
                    },
                }

        


        self.SysDict_rateshift = {
                "lumi_13TeV_2016":{
                            "process": ["Signal", "TT","sT","TT_untagged","sT_untagged", "ttV","VV", "DY"],    
                            "description":"luminosity uncertainty",
                            "isShape":False,
                            "channel":["MuMu","MuEl","ElEl"],
                            "value": 1.025,
                        },
                "ttbar_xsec":{
                        "process": ["TT", "TT_untagged"],
                        "description":"TTbar cross section uncertainty",
                        "channel": ["MuMu","MuEl","ElEl"],
                        "isShape":False,
                        "value": 1.053,
                        },
                "singleTop_xsec":{
                        "process":["sT", "sT_untagged"],
                        "description":"singletop cross section uncertainty",
                        "channel":["MuMu","MuEl","ElEl"],
                        "isShape":False,
                        "value":1.072,
                        },
                "dy_mc_xsec":{
                        "process":["DY"],
                        "description":"Drell_Yan cross section uncertainty",
                        "channel":["MuEl"],
                        "isShape":False,
                        "value":1.05,
                        },
                "dy_rwgt_norm_MuMu":{
                        "process":["TT_untagged","sT_untagged","data_untagged"],
                        "channel":["MuMu"],
                        "description":"Data-driven estimation normalization uncertainty",
                        "isShape":False,
                        "value":1.05,
                        },
                "dy_rwgt_norm_ElEl":{
                        "process":["TT_untagged","sT_untagged","data_untagged"],
                        "channel":["ElEl"],
                        "description":"Data-driven estimation normalization uncertainty",
                        "isShape":False,
                        "value":1.05,
                        },
                }


        self.addRateShift = True

        self.filedict = inputfiledict

        self.nominalWeight = nominalWeight

        self.systematiclist = systematiclist

        self.DYdatadriven = DYdatadriven
        self.shortnames_datadriven = ["TT_untagged", "sT_untagged"]

        

        self.treename = "evtreeHME_nn"
        #self.treename_signal = "Friends"
        self.treename_signal = "evtreeHME_nn"

    #def initialize1D(self, inputfiledict, nominalWeight, systematiclist, todraw, xtitle, xbins, cuts, DYdatadriven, untagged_filedict):
    def initialize1D(self, todraw, xtitle, xbins, cuts):
        """Setup some initial variables"""
        print "============    initialize1D  =================================="
        print "1D shape systematics ", todraw, " x ",xtitle, " cuts ",cuts
        self.todraw = todraw
        self.xtitle = xtitle
        self.xbins = xbins
        self.cuts = cuts
        #self.shortname = shortname 

        self.sample_systematic_hist = {}
        self.systematic_hist = {}
        self.finalhist = {}
        self.channel_shortname_systematic_hist = {}
        for channel in self.channels:
            self.channel_shortname_systematic_hist[channel] = {}

        print "============ end of   initialize1D  =================================="

    def initialize2D(self, todraw, xtitle, xbins, ytitle, ybins, cuts):
        """Setup some initial variables"""
        print "============    initialize2D  =================================="
        print "2D shape systematics ", todraw, " x ",xtitle, " y ",ytitle," cuts ",cuts
        self.todraw2D = todraw
        self.xtitle = xtitle
        self.xbins = xbins
        self.cuts = cuts
        self.ytitle = ytitle
        self.ybins = ybins
        #self.shortname = shortname 

        self.sample_systematic_hist = {}
        self.systematic_hist = {}
        self.finalhist = {}
        self.channel_shortname_systematic_hist = {}
        for channel in self.channels:
            self.channel_shortname_systematic_hist[channel] = {}

        print "============ end of   initialize2D  =================================="

    def get_xsection_eventweightsum_file(self, filename):
        tfile = ROOT.TFile(filename, "READ")
        xsec = tfile.Get("cross_section")
        event_weight_sum = tfile.Get("event_weight_sum")
        return xsec.GetVal(),event_weight_sum.GetVal()

    def get_xsection_eventweightsum_tree(self, filename, treename):
        tree = ROOT.TChain( treename )
        tree.Add(filename)
        n = tree.GetEntry()
        i = 0
        xsec = 0.0; event_weight_sum =0.0;
        for i in range(0, 100):
            tree.GetEntry(i)
            cross_section = tree.cross_section
            weisum = tree.event_weight_sum
            if i == 0:
                xsec = cross_section
                event_weight_sum = weisum
            else:
                if abs(xsec-cross_section)>0.01*xsec or abs(event_weight_sum - weisum)> 0.01*event_weight_sum:
                    print "WARNING: cross_section or event_weight_sum may be not a single value, xsec  ", xsec," event_weight_sum ",event_weight_sum," this entry ",cross_section," ",weisum 
        return xsec,event_weight_sum

    def getTotalRateShift(self, shortname, channel):
        delta_shift_square = 0.0
        for sys in self.SysDict_rateshift:
            if channel in self.SysDict_rateshift[sys]["channel"] and (shortname in self.SysDict_rateshift[sys]["process"] or ("Radion" in shortname and "Signal" in self.SysDict_rateshift[sys]["process"])):
                delta = self.SysDict_rateshift[sys]["value"] -1.0
                delta_shift_square = delta_shift_square + delta * delta
        print "Sample ",shortname, " channel ",channel, " total rateshift ", 1.0 + sqrt(delta_shift_square)
        return 1.0 + sqrt(delta_shift_square)


    def GetQCDScaleHist(self, filename):
        tfile = ROOT.TFile(filename, "READ")
        hist = tfile.Get("CountWeightedLHEWeightScale")
        finalhist = hist.Clone()
        finalhist.SetDirectory(0)
        return finalhist


    #def QCDScaleSystematic(self, chain, weight, nominal_event_weight_sum, QCDscalehist, plotname):
    def QCDScaleSystematic(self, chain, weight, plotname):

        
        #####
        ####LHE scale variation weights (w_var / w_nominal); [0] is renscfact=0.5d0 facscfact=0.5d0 ; [1] is renscfact=0.5d0 facscfact=1d0 ; [2] is renscfact=0.5d0 facscfact=2d0 ; [3] is renscfact=1d0 facscfact=0.5d0 ; [4] is renscfact=1d0 facscfact=1d0 ; [5] is renscfact=1d0 facscfact=2d0 ; [6] is renscfact=2d0 facscfact=0.5d0 ; [7] is renscfact=2d0 facscfact=1d0 ; [8] is renscfact=2d0 facscfact=2d0
        ####
        indexlist = [0, 1, 3, 5, 7, 8]
        QCDscale_allshapes = []
        for index in indexlist:
            #this_event_weight_sum = QCDscalehist.GetBinContent(index + 1)
            #finalcut = "LHEScaleWeight[%d]"%index +"*"+weight 
            finalcut = "LHEScaleWeight_%d"%index +"*"+weight 
            hist = ROOT.TH1F("QCDScale_shape%d"%index, "QCDScale_shape%d"%index, len(self.xbins)-1, self.xbins)
            QCDscale_allshapes.append(hist)
            chain.Draw(self.todraw + ">> QCDScale_shape%d"%index, finalcut)
            #print "finalcut ",finalcut, " hist integral ", hist.Integral()

        hist_up = ROOT.TH1F("QCDscaleup","QCDscaleup", len(self.xbins)-1, self.xbins)
        hist_down = ROOT.TH1F("QCDscaledown", "QCDscaledown", len(self.xbins)-1, self.xbins)
        for bin in xrange( hist_up.GetNbinsX()):
            binvalues = []
            for shape in QCDscale_allshapes:
                binvalues.append(shape.GetBinContent(bin + 1))
            hist_up.SetBinContent(bin+1, max(binvalues))
            hist_down.SetBinContent(bin+1, min(binvalues))
            #print "bin ",bin," binvalues ",binvalues, " max ", max(binvalues), " min ",min(binvalues)


        #print "final integral up ", hist_up.Integral(), " down ", hist_down.Integral(), " nominal ",QCDscale_allshapes[3].Integral()
        plotQCDScale = False
        if plotQCDScale:
            samplename = plotname.split('/')[-1]
            hs = ROOT.THStack(samplename, "QCD scale uncertainty, intermedian plot")
            #hist_up.SetMarkerColor(2)
            #hist_up.SetMarkerStyle(22)
            #hist_down.SetMarkerColor(4)
            #hist_down.SetMarkerStyle(23)
            hist_up.SetLineColor(1)
            #hist_up.SetLineStyle(2)
            hist_up.SetLineWidth(2)
            hist_down.SetLineColor(1)
            #hist_down.SetLineStyle(2)
            hist_down.SetLineWidth(2)
            colors = [800-4, 820-3, 900-3, 860-3, 616-7, 432+2, 400+2]
            legend = ROOT.TLegend(0.74,0.5,0.84,0.5+7*.05); 
            legend.SetTextSize(0.04); legend.SetTextFont(42)
            for i, shape in enumerate(QCDscale_allshapes):
                shape.SetLineColor(colors[i])
                hs.Add(shape)
                legend.AddEntry(shape, "QCDscale_%d"%indexlist[i],"l")
            legend.SetBorderSize(0)
            legend.AddEntry(hist_up, "QCDscale Up","l")
            legend.AddEntry(hist_down, "QCDscale down","l")
     
            scale_c = ROOT.TCanvas("scale", "scale",800, 600)
            hs.Draw("hist nostack")
            hist_up.Draw("histsame")
            hist_down.Draw("histsame")
            legend.Draw("same")
            hs.GetHistogram().GetXaxis().SetTitle(self.xtitle)
            hs.GetHistogram().GetYaxis().SetTitle("Events")
            tex1 = ROOT.TLatex(0.13,0.87, samplename)
            tex1.SetNDC(); tex1.SetTextSize(.035)
            tex1.Draw("same")
            #plotdir = "DataDriven_DY_plots/"
            scale_c.SaveAs(plotname+"_QCDscale.pdf")
          


        hist_up.SetDirectory(0)
        hist_down.SetDirectory(0)
        return hist_up,hist_down



    #################################################################################################
    #### do systematics for one file
    #################################################################################################
    def plotSystematics_singleFile(self, key, channel, samplename, filename, nominalWeight): 

        
        self.sample_systematic_hist[samplename] = {}

        #xsec,event_weight_sum = self.get_xsection_eventweightsum_file(filename)
        treename = self.treename 
        if "Radion" in samplename :
            treename = self.treename_signal
        xsec,event_weight_sum = self.get_xsection_eventweightsum_tree(filename, treename)
        weight = nominalWeight+"*{cross_section}*1000.0/{event_weight_sum}".format(cross_section = xsec, event_weight_sum = event_weight_sum)
        reweight = xsec*1000.0/event_weight_sum
        #print "sample ",samplename, " xsec ",xsec, " event_weight_sum ",event_weight_sum, " reweight ",reweight
        #weight = nominalWeight+"*cross_section*1000.0/event_weight_sum"
        #if "Radion" in key or "Graviton" in key:
        #    #xsec = 5.0 #pb
        #    xsec = 1.0e-3 # pb = 1fb
        chain = ROOT.TChain(self.treename)
        if "Radion" in samplename:
            chain = ROOT.TChain(self.treename_signal)
        chain.Add( filename )
        self.sample_systematic_hist[samplename]["nominal"] = ROOT.TH1F(key+"_"+channel+"_%s"%samplename+"_nominal", key+"_"+channel+"_%s"%samplename+"_nominal", len(self.xbins)-1, self.xbins)
        self.sample_systematic_hist[samplename]["nominal_noweight"] = ROOT.TH1F(key+"_"+channel+"_%s"%samplename+"_nominal_noweight", key+"_"+channel+"_%s"%samplename+"_nominal_noweight", len(self.xbins)-1, self.xbins) ## not possible
        #self.sample_systematic_hist[samplename]["weight"] = ROOT.TH1F(key+"_"+channel+"_%s"%samplename+"_weight", key+"_"+channel+"_%s"%samplename+"_weight", len(self.xbins)-1, self.xbins)
        finalcut = "("+ self.cuts +"&&"+ self.channelcuts[channel] +")*"+weight
        #print "sample ",samplename, " xsec ",xsec, "  event_weight_sum ",event_weight_sum, " weight ",weight , " finalcut ",finalcut
        chain.Draw(self.todraw + ">> " + self.sample_systematic_hist[samplename]["nominal"].GetName(), finalcut)
        chain.Draw(self.todraw + ">> " + self.sample_systematic_hist[samplename]["nominal_noweight"].GetName(), self.cuts +"&&"+ self.channelcuts[channel])
        #chain.Draw(  "1 >> " + self.sample_systematic_hist[samplename]["weight"].GetName(), finalcut)
        #print "nominal integral ", self.sample_systematic_hist[samplename]["nominal"].Integral()
        self.sample_systematic_hist[samplename]["nominal"].SetDirectory(0)
        self.sample_systematic_hist[samplename]["nominal_noweight"].SetDirectory(0)
        #self.sample_systematic_hist[samplename]["weight"].SetDirectory(0)


        ###now plotting systematics
        for sys in self.systematiclist:
            self.sample_systematic_hist[samplename][sys] = {}

            if sys == "QCDscale":
                faileddatasets =  ["WWToLNuQQ_aTGC_13TeV-madgraph-pythia8", "ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1"]
                faileddatasets.append("ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1")
                if (samplename in faileddatasets) or (self.DYdatadriven and samplename.replace("_untagged","") in faileddatasets):
                    self.sample_systematic_hist[samplename][sys]["up"] = None
                    self.sample_systematic_hist[samplename][sys]["down"] = None
                    continue
                plotdir = "QCDScale_intermedianplots/"
                #QCDscalehist = self.GetQCDScaleHist(filename)
                plotname = os.path.join(plotdir, samplename+"_"+channel)
                weight_cuts = "("+ self.cuts +"&&"+ self.channelcuts[channel] +")*"+weight
                #hist_up, hist_down = self.QCDScaleSystematic(chain, weight_cuts, event_weight_sum, QCDscalehist, plotname)
                hist_up, hist_down = self.QCDScaleSystematic(chain, weight_cuts, plotname)
                hist_up.SetName(key+"_"+channel+"_%s"%samplename+"_QCDscaleup")
                hist_down.SetName(key+"_"+channel+"_%s"%samplename+"_QCDscaledown")
                #print "samplename ",samplename, " up ",hist_up.Print("ALL")," down ",hist_down.Print("ALL")," nominal ",self.sample_systematic_hist[samplename]["nominal"].Print("ALL")
                self.sample_systematic_hist[samplename][sys]["up"] = hist_up
                self.sample_systematic_hist[samplename][sys]["down"] = hist_down
                continue

            if sys == "MC_statistical":
                hist_up = ROOT.TH1F(key+"_"+channel+"_%s"%samplename+"_"+sys+"_up", key+"_"+channel+"_%s"%samplename+"_"+sys+"_up", len(self.xbins)-1, self.xbins)
                hist_down = ROOT.TH1F(key+"_"+channel+"_%s"%samplename+"_"+sys+"_down", key+"_"+channel+"_%s"%samplename+"_"+sys+"_down", len(self.xbins)-1, self.xbins)
                ### why it is 2D here?
                #print("workout the systematic uncertainty for histogram ",self.sample_systematic_hist[samplename]["nominal"].GetName()," GetNbinsX ", self.sample_systematic_hist[samplename]["nominal"].GetNbinsX(), " GetNbinsY ", self.sample_systematic_hist[samplename]["nominal"].GetNbinsY())
                for bin in xrange(self.sample_systematic_hist[samplename]["nominal"].GetNbinsX()):
                    #for ybin in xrange(self.sample_systematic_hist[samplename]["nominal"].GetNbinsY()):
                    binvalue = self.sample_systematic_hist[samplename]["nominal"].GetBinContent(bin+1)
                    rawbinvalue = self.sample_systematic_hist[samplename]["nominal_noweight"].GetBinContent(bin+1)
                    weight_thisbin = 1.0
                    if rawbinvalue > 0:
                        weight_thisbin = binvalue/rawbinvalue
                    err = sqrt(abs(rawbinvalue))*weight_thisbin
                    hist_up.SetBinContent(bin+1, binvalue + err)
                    hist_down.SetBinContent(bin+1,  binvalue - err)
                    ###self.sample_systematic_hist[samplename]["nominal"].SetBinError(bin+1, ybin+1, err)
                    if err >  0.1*abs(binvalue) and (key == "TT"):
                        print "warning, in this bin ",bin+1," stat_err ",err, " binvalue ",binvalue, " is larger than 10% of binvalue in TTbar!!!"
                    #if bin == 1:
                    #    print(key+"_"+channel+"_%s"%samplename," binvalue ", binvalue," binvalue_raw ", rawbinvalue," bin_weighthist ", weight_thisbin, " staterr ", err)

                self.sample_systematic_hist[samplename][sys]["up"] = hist_up
                self.sample_systematic_hist[samplename][sys]["down"] = hist_down
                continue

            for plottype in ["up","down"]:

                suffix = sys+"_"+plottype
                hist = ROOT.TH1F(key+"_"+channel+"_%s"%samplename+"_%s"%(suffix), key+"_"+channel+"_%s"%samplename+"_%s"%(suffix), len(self.xbins)-1, self.xbins)
                self.sample_systematic_hist[samplename][sys][plottype] = hist
                
                sysweight = self.SysDict[sys]["weight"].replace("nominal", plottype)
                thisweight = sysweight+"*"+weight
                finalcut = "(" + self.cuts +"&&"+ self.channelcuts[channel] + ")*"+thisweight
                #chlist[samplename].Draw(todraw + ">> " + self.sample_systematic_hist[samplename][sys][plottype].GetName(), finalcut)
                #print "sample ",samplename," file ",self.filedict[key][samplename]['path']," weight ",thisweight," todraw ",self.todraw," finalcut ",finalcut
                chain.Draw(self.todraw + ">> " + self.sample_systematic_hist[samplename][sys][plottype].GetName(), finalcut)
                self.sample_systematic_hist[samplename][sys][plottype].SetDirectory(0)
                

    #################################################################################################
    ##plot systematcs of all sub-processes in one process like TT, sT, DY...
    #################################################################################################
    def plotSystematics(self, key, channel):
        
        self.shortname = key
        self.hist_data_untagged = None
        self.sample_systematic_hist = {}
       
        #print "plotSystematics ",key, " channel ", channel
        if self.shortname == "DY" and self.DYdatadriven and (channel == "MuMu" or channel == "ElEl"):
            ##plot data
            self.hist_data_untagged = ROOT.TH1F("data_untagged_"+channel, "data_untagged_"+channel, len(self.xbins)-1, self.xbins)
            Mbtag_weight = "dy_Mbtag_weight"
            untagged_suffix = "_untagged"
            data_untagged_name = "DoubleMuon"
            if channel == "ElEl":
                data_untagged_name = "DoubleEG"
            datafile_untagged =  self.filedict["Data"+untagged_suffix][data_untagged_name]['path']
            ch_d = ROOT.TChain(self.treename)
            ch_d.AddFile(datafile_untagged)
            cut_data = "("+ self.cuts +"&&"+ self.channelcuts[channel] +")*"+Mbtag_weight
            ch_d.Draw(self.todraw + ">> " + self.hist_data_untagged.GetName(), cut_data)
            self.hist_data_untagged.SetDirectory(0)

            weight_datadriven = self.nominalWeight+"* "+ Mbtag_weight
            for key_datadriven in self.shortnames_datadriven:
                for samplename_datadriven in self.filedict[key_datadriven].keys():
                    self.plotSystematics_singleFile(key, channel, samplename_datadriven, self.filedict[key_datadriven][samplename_datadriven]['path'], weight_datadriven)
        else:
            for iname, samplename in enumerate(self.filedict[self.shortname].keys()):
                self.plotSystematics_singleFile(self.shortname, channel, samplename, self.filedict[self.shortname][samplename]['path'], self.nominalWeight)

    
    #################################################################################################
    ##combine subprocess systematcs for one process: like TT, sT, ttV, TT_untagged, sT_untagged
    #################################################################################################
    def combineSystematics(self, shortname, channel):
        ##combine to final plots
        self.shortname = shortname
        suffix = "nominal"
        #print "combineSystematics ", shortname, " channel ", channel
        self.finalhist["nominal"] = ROOT.TH1F(self.shortname+"_"+channel, self.shortname+"_"+channel, len(self.xbins)-1, self.xbins)
        allsamples_combine = None
        if "untagged" in self.shortname and self.DYdatadriven and (channel == "MuMu" or channel == "ElEl"):
            allsamples_combine = self.filedict[shortname]
        else:
            allsamples_combine = self.filedict[self.shortname].keys()

        for iname, samplename in enumerate(allsamples_combine):
            self.finalhist["nominal"].Add( self.sample_systematic_hist[samplename]["nominal"] )
        self.finalhist["nominal"].SetDirectory(0)

        rateshift = 1.0
        if self.addRateShift:
            rateshift = self.getTotalRateShift(shortname, channel)

        dummyplots = {}
        totalSys = {}
        for sys in self.systematiclist:
            self.systematic_hist[sys] = {}
            totalSys[sys] = {}
        for plottype in ["up","down"]:

            dummyplots[plottype] = self.cloneDummyHistogram( self.finalhist["nominal"] )
            dummyplots[plottype].SetName("dummy"+plottype)
            for sys in self.systematiclist:
                totalSys[sys][plottype] = 0.0
                suffix = "_"+sys+"_"+plottype
                #hist = ROOT.TH1F(self.shortname+"_"+channel+suffix, self.shortname+"_"+channel+suffix, len(self.xbins)-1, self.xbins)
                self.systematic_hist[sys][plottype] = self.cloneDummyHistogram( self.finalhist["nominal"] )
                self.systematic_hist[sys][plottype].SetName(self.shortname+"_"+channel+suffix)
                #print "systematic hist name ",self.shortname+"_"+channel+suffix
                for iname, samplename in enumerate(allsamples_combine):
                    if self.sample_systematic_hist[samplename][sys][plottype]:
                        self.systematic_hist[sys][plottype].Add( self.sample_systematic_hist[samplename][sys][plottype] )
                    else:
                        ##if sys shift not available, still need to add the nominal 
                        self.systematic_hist[sys][plottype].Add( self.sample_systematic_hist[samplename]["nominal"] )
                #print "plottype ",plottype, self.systematic_hist[sys][plottype].Print("ALL")," nominal ",self.finalhist["nominal"].Print("ALL")
                self.systematic_hist[sys][plottype].SetDirectory(0)
                for bin in xrange(self.systematic_hist[sys][plottype].GetNbinsX()):
                    sysvalue = self.systematic_hist[sys][plottype].GetBinContent(bin+1) - self.finalhist["nominal"].GetBinContent(bin+1) ##error: hist_up - hist_nominal
                    #print "plottype ",plottype, " bin ",bin," systematic_hist ",self.systematic_hist[sys][plottype].GetBinContent(bin+1)," nominal ",self.finalhist["nominal"].GetBinContent(bin+1)  ," sysvalue ",sysvalue ," dummyplot bin content ", dummyplots[plottype].GetBinContent(bin+1)
                    totalSys[sys][plottype] = totalSys[sys][plottype] + sysvalue
                    bincontent = dummyplots[plottype].GetBinContent(bin+1) + sysvalue * sysvalue
                    dummyplots[plottype].SetBinContent(bin+1, bincontent)
                totalSys[sys][plottype] = abs(totalSys[sys][plottype])/self.finalhist["nominal"].Integral()

            for bin in xrange( dummyplots[plottype].GetNbinsX()):
                bincontent = dummyplots[plottype].GetBinContent(bin+1) 
                nominal =  self.finalhist["nominal"].GetBinContent(bin+1)
                finalError2 = bincontent + nominal*nominal*(rateshift-1.0)*(rateshift-1.0)
                dummyplots[plottype].SetBinContent(bin+1, sqrt(finalError2))


            ### change the histname and change the statistical uncertainty!!!, for 1D and 2D, Tao
            self.finalhist[plottype] = self.finalhist["nominal"].Clone()
            self.finalhist[plottype].SetName(self.shortname+"_"+channel+plottype)
            addweight = 1.0
            #print "plottype ", plottype, " dummpyplot ",dummyplots[plottype].Print("ALL")
            if plottype == "down":
                addweight = -1.0
            self.finalhist[plottype].Add(dummyplots[plottype], addweight)
            self.finalhist[plottype].SetDirectory(0)

        #### end of plot type loop

        #print "nominal ", self.finalhist["nominal"].Print("ALL")
        ##add one side 
        self.finalhist["one_sided"] = self.cloneDummyHistogram( self.finalhist["nominal"] )
        self.finalhist["one_sided"].SetName(self.shortname+"_"+channel+"one_sided")
        self.finalhist["nominal_allSys"] = self.cloneDummyHistogram( self.finalhist["nominal"] )
        self.finalhist["nominal_allSys"].SetName(self.shortname+"_"+channel+"allSys")
        
        for bin in xrange(self.finalhist["one_sided"].GetNbinsX()):
            err_up = self.finalhist["up"].GetBinContent(bin + 1)
            err_down = self.finalhist["down"].GetBinContent(bin + 1)
            staterr_up = self.systematic_hist["MC_statistical"]["up"].GetBinContent(bin+1)
            nominal = self.finalhist["nominal"].GetBinContent(bin + 1)
            onesided = (abs(nominal - err_up) + abs(nominal - err_down))/2.0
            staterr = abs(staterr_up-nominal) ## should be symetric
            #print "err_up ",err_up, " err_down ",err_down, " nominal ",nominal," oneside ",onesided," error/nominal ", onesided/nominal
            #sys_stat_total = sqrt(onesided*onesided + nominal)
            ## final nominal shape with what uncertainty ??
            ###Attention!!
            self.finalhist["one_sided"].SetBinContent(bin+1, nominal + onesided)
            self.finalhist["nominal_allSys"].SetBinContent(bin+1, nominal) 
            self.finalhist["nominal_allSys"].SetBinError(bin+1, onesided) ## including all uncertainties
            self.finalhist["nominal"].SetBinError(bin+1, staterr) ## only including stat uncertainties
            ####debug
            #if bin == 1:
            #    print(self.shortname+"_"+channel, " bin1 ", nominal, " MC_statistical_up ",staterr_up," staterr ", staterr, " allSys ",  onesided)


        self.systematic_hist["nominal"] = self.finalhist["nominal"].Clone()
        #print "Integral nominal ", self.finalhist["nominal"].Integral(), " up ",self.finalhist["up"].Integral()," down ",self.finalhist["down"].Integral(), " one sided ",self.finalhist["one_sided"].Integral()
        #print "totalSys ",totalSys," self.systematic_hist ",self.systematic_hist
        for sys in self.systematiclist:
            totalSys[sys]["oneside"] = (totalSys[sys]["up"] + totalSys[sys]["down"])/2.0
            print "Systematic ",sys," samplename ",self.shortname, " up ",totalSys[sys]["up"], " down " ,totalSys[sys]["down"]," oneside ",totalSys[sys]["oneside"]
        print "Final ", self.finalhist["nominal"].Integral(), " one_sided ", self.finalhist["one_sided"].Integral()," error ", self.finalhist["one_sided"].Integral()/self.finalhist["nominal"].Integral()-1.0


    #################################################################################################
    ### 
    #################################################################################################
    def runSystematics(self, shortname, channel):
        ### step 1, get histogram from single file
        ### step 2, get histogram for one type of sample: TTbar, ST, VV, ttV , DY
        self.plotSystematics(shortname, channel)
        if shortname == "DY" and self.DYdatadriven and (channel == "MuMu" or channel == "ElEl"):
            datadriven_finalhist ={}
            for name in self.shortnames_datadriven:
                self.combineSystematics(name, channel)
                datadriven_finalhist[name] = self.finalhist.copy()
                self.channel_shortname_systematic_hist[channel][name] = self.systematic_hist.copy()
                #print "Data-driven ", datadriven_finalhist[name]["nominal"].Print("ALL")

            ### calculate final hist for data-driven DY    
            self.finalhist = {}
            for plottype in ["nominal","up","down","one_sided","nominal_allSys"]:
                suffix = "_"+plottype
                if plottype == "nominal":
                    suffix = ""
                self.finalhist[plottype] = self.cloneDummyHistogram(self.hist_data_untagged)
                self.finalhist[plottype].SetName(shortname + "_"+ channel+ suffix)

            for  bin in xrange(self.hist_data_untagged.GetNbinsX()):
                err_up_square = 0.0
                err_down_square = 0.0
                nominal_mc = 0.0
                data_binvalue = self.hist_data_untagged.GetBinContent(bin+1)
                for name in self.shortnames_datadriven:
                    thisnominal = datadriven_finalhist[name]["nominal"].GetBinContent(bin+1)
                    thisup = datadriven_finalhist[name]["up"].GetBinContent(bin+1)
                    thisdown = datadriven_finalhist[name]["down"].GetBinContent(bin+1)
                    err_up_square = err_up_square + (thisnominal-thisup)*(thisnominal-thisup)
                    err_down_square = err_down_square + (thisnominal-thisdown)*(thisnominal-thisdown)
                    nominal_mc  = nominal_mc+thisnominal
                one_sided_err = (sqrt(err_up_square)+sqrt(err_down_square))/2.0
                self.finalhist["nominal"].SetBinContent(bin+1, data_binvalue - nominal_mc)
                self.finalhist["nominal_allSys"].SetBinContent(bin+1, data_binvalue - nominal_mc)
                self.finalhist["up"].SetBinContent(bin+1, data_binvalue - nominal_mc + sqrt(err_up_square))
                self.finalhist["down"].SetBinContent(bin+1, data_binvalue - nominal_mc - sqrt(err_down_square))
                self.finalhist["one_sided"].SetBinContent(bin+1, data_binvalue - nominal_mc + one_sided_err)
                self.finalhist["nominal_allSys"].SetBinError(bin+1, one_sided_err)
                ## what is stat err for data -driven??
            #print "Data-Driven DY final ", self.finalhist["nominal"].Print("ALL")


        else:
            self.combineSystematics(shortname, channel)
            self.channel_shortname_systematic_hist[channel][shortname] = self.systematic_hist.copy()
        #print "shortname ",shortname, " all process systematics ", self.channel_shortname_systematic_hist

    def addTH1withError(self, hist1, hist2, c2=1.0):
        
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
            
    def cloneDummyHistogram(self,rHistogram):
        """Clone a histogram and set content/error to 0 -- for 1-sided systematics"""
        h_dummy = rHistogram.Clone()
        h_dummy.SetName("dummy")
        for bin in xrange(h_dummy.GetNbinsX()):
            h_dummy.SetBinContent(bin+1,0.0)
            h_dummy.SetBinError(bin+1,0.0)

        return h_dummy


    def cloneDummyHistogram2D(self,rHistogram):
        """Clone a histogram and set content/error to 0 -- for 1-sided systematics"""
        h_dummy = rHistogram.Clone()
        h_dummy.SetName("dummy")
        for bin in xrange(h_dummy.GetNbinsX()):
            for biny in xrange(h_dummy.GetNbinsY()):
                h_dummy.SetBinContent(bin+1, biny+1,0.0)
                h_dummy.SetBinError(bin+1, biny+1, 0.0)

        return h_dummy



    def writeSystematicsToFile(self, directory):

        ## how to write histogram into file 
        self.outfile = os.path.join(directory, self.shortname+"_systematics.root")
        print "self.outfile ", self.outfile
        tfile = ROOT.TFile(self.outfile, "recreate")
        for samplename in self.sample_systematic_hist.keys():
            #self.sample_systematic_hist[samplename]["nominal"].Print("ALL")
            self.sample_systematic_hist[samplename]["nominal"].SetDirectory(tfile)
            self.sample_systematic_hist[samplename]["nominal"].Write()
            for sys in self.sample_systematic_hist[samplename].keys():
                for plottype in ["up","down"]:
                    self.sample_systematic_hist[samplename][sys][plottype].SetDirectory(tfile)
                    self.sample_systematic_hist[samplename][sys][plottype].Write()

        self.systematic_hist["nominal"].SetDirectory(tfile)
        self.systematic_hist["nominal"].Write()
        for sys in self.systematic_hist.keys():
            for plottype in ["up","down"]:
                self.systematic_hist[sys][plottype].SetDirectory(tfile)
                self.systematic_hist[sys][plottype].Write()

        """
        for plotytype in ["nominal", "up","down"]
            self.finalhist[plottype].SetDirectory(tfile)
            self.finalhist[plottype].Write()
        """

        tfile.Close()







    #################################################################################################
    ###  2D systematics starts from  here
    #################################################################################################

    #### QCD scale for 2D
    def QCDScaleSystematic2D(self, chain, weight, plotname):

        
        #####
        ####LHE scale variation weights (w_var / w_nominal); [0] is renscfact=0.5d0 facscfact=0.5d0 ; [1] is renscfact=0.5d0 facscfact=1d0 ; [2] is renscfact=0.5d0 facscfact=2d0 ; [3] is renscfact=1d0 facscfact=0.5d0 ; [4] is renscfact=1d0 facscfact=1d0 ; [5] is renscfact=1d0 facscfact=2d0 ; [6] is renscfact=2d0 facscfact=0.5d0 ; [7] is renscfact=2d0 facscfact=1d0 ; [8] is renscfact=2d0 facscfact=2d0
        ####
        indexlist = [0, 1, 3, 5, 7, 8]
        QCDscale_allshapes = []
        for index in indexlist:
            #this_event_weight_sum = QCDscalehist.GetBinContent(index + 1)
            #finalcut = "LHEScaleWeight[%d]"%index +"*"+weight 
            finalcut = "LHEScaleWeight_%d"%index +"*"+weight 
            hist = ROOT.TH2F("QCDScale_shape%d"%index, "QCDScale_shape%d"%index, len(self.xbins)-1, self.xbins, len(self.ybins) -1, self.ybins)
            QCDscale_allshapes.append(hist)
            chain.Draw(self.todraw2D + ">> QCDScale_shape%d"%index, finalcut)
            #print "finalcut ",finalcut, " hist integral ", hist.Integral()

        hist_up = ROOT.TH2F("QCDscaleup","QCDscaleup", len(self.xbins)-1, self.xbins,  len(self.ybins) -1, self.ybins)
        hist_down = ROOT.TH2F("QCDscaledown", "QCDscaledown", len(self.xbins)-1, self.xbins,  len(self.ybins) -1, self.ybins)
        for bin in xrange( hist_up.GetNbinsX()):
            for ybin in range( hist_up.GetNbinsY()):
                binvalues = []
                for shape in QCDscale_allshapes:
                    binvalues.append(shape.GetBinContent(bin + 1, ybin+1))
                hist_up.SetBinContent(bin+1, ybin+1, max(binvalues))
                hist_down.SetBinContent(bin+1, ybin+1, min(binvalues))
            #print "bin ",bin," binvalues ",binvalues, " max ", max(binvalues), " min ",min(binvalues)


        #print "final integral up ", hist_up.Integral(), " down ", hist_down.Integral(), " nominal ",QCDscale_allshapes[3].Integral()
        #plotQCDScale = False
        #if plotQCDScale:
        #    samplename = plotname.split('/')[-1]
        #    hs = ROOT.THStack(samplename, "QCD scale uncertainty, intermedian plot")
        #    #hist_up.SetMarkerColor(2)
        #    #hist_up.SetMarkerStyle(22)
        #    #hist_down.SetMarkerColor(4)
        #    #hist_down.SetMarkerStyle(23)
        #    hist_up.SetLineColor(1)
        #    #hist_up.SetLineStyle(2)
        #    hist_up.SetLineWidth(2)
        #    hist_down.SetLineColor(1)
        #    #hist_down.SetLineStyle(2)
        #    hist_down.SetLineWidth(2)
        #    colors = [800-4, 820-3, 900-3, 860-3, 616-7, 432+2, 400+2]
        #    legend = ROOT.TLegend(0.74,0.5,0.84,0.5+7*.05); 
        #    legend.SetTextSize(0.04); legend.SetTextFont(42)
        #    for i, shape in enumerate(QCDscale_allshapes):
        #        shape.SetLineColor(colors[i])
        #        hs.Add(shape)
        #        legend.AddEntry(shape, "QCDscale_%d"%indexlist[i],"l")
        #    legend.SetBorderSize(0)
        #    legend.AddEntry(hist_up, "QCDscale Up","l")
        #    legend.AddEntry(hist_down, "QCDscale down","l")
     
        #    scale_c = ROOT.TCanvas("scale", "scale",800, 600)
        #    hs.Draw("hist nostack")
        #    hist_up.Draw("histsame")
        #    hist_down.Draw("histsame")
        #    legend.Draw("same")
        #    hs.GetHistogram().GetXaxis().SetTitle(self.xtitle)
        #    hs.GetHistogram().GetYaxis().SetTitle("Events")
        #    tex1 = ROOT.TLatex(0.13,0.87, samplename)
        #    tex1.SetNDC(); tex1.SetTextSize(.035)
        #    tex1.Draw("same")
        #    #plotdir = "DataDriven_DY_plots/"
        #    scale_c.SaveAs(plotname+"_QCDscale.pdf")
          


        hist_up.SetDirectory(0)
        hist_down.SetDirectory(0)
        return hist_up,hist_down



    #################################################################################################
    #### do systematics for one file, 2D
    #################################################################################################
    def plotSystematics_singleFile2D(self, key, channel, samplename, filename, nominalWeight): 

        
        self.sample_systematic_hist[samplename] = {}

        #xsec,event_weight_sum = self.get_xsection_eventweightsum_file(filename)
        treename = self.treename 
        if "Radion" in samplename:
            treename = self.treename_signal
        xsec,event_weight_sum = self.get_xsection_eventweightsum_tree(filename, treename)
        weight = nominalWeight+"*{cross_section}*1000.0/{event_weight_sum}".format(cross_section = xsec, event_weight_sum = event_weight_sum)
        reweight = xsec*1000.0/event_weight_sum
        #print "sample ",samplename, " xsec ",xsec, " event_weight_sum ",event_weight_sum, " reweight ",reweight," filename ",filename
        #weight = nominalWeight+"*cross_section*1000.0/event_weight_sum"
        #if "Radion" in key or "Graviton" in key:
        #    #xsec = 5.0 #pb
        #    xsec = 1.0e-3 # pb = 1fb
        chain = ROOT.TChain(self.treename)
        if "Radion" in samplename:
            chain = ROOT.TChain(self.treename_signal)
        chain.Add( filename )
        self.sample_systematic_hist[samplename]["nominal"] = ROOT.TH2F(key+"_"+channel+"_%s"%samplename+"_nominal", key+"_"+channel+"_%s"%samplename+"_nominal", len(self.xbins)-1, self.xbins,  len(self.ybins) -1, self.ybins)
        self.sample_systematic_hist[samplename]["nominal_noweight"] = ROOT.TH2F(key+"_"+channel+"_%s"%samplename+"_nominal_noweight", key+"_"+channel+"_%s"%samplename+"_nominal_noweight", len(self.xbins)-1, self.xbins,  len(self.ybins) -1, self.ybins)
        finalcut = "("+ self.cuts +"&&"+ self.channelcuts[channel] +")*"+weight
        #print "todraw ",self.todraw2D," sample ",samplename, " xsec ",xsec, "  event_weight_sum ",event_weight_sum, " weight ",weight , " finalcut ",finalcut
        chain.Draw(self.todraw2D + ">> " + self.sample_systematic_hist[samplename]["nominal"].GetName(), finalcut)
        chain.Draw(self.todraw2D + ">> " + self.sample_systematic_hist[samplename]["nominal_noweight"].GetName(), self.cuts +"&&"+ self.channelcuts[channel])
        #print "nominal integral ", self.sample_systematic_hist[samplename]["nominal"].Integral()
        self.sample_systematic_hist[samplename]["nominal"].SetDirectory(0)
        self.sample_systematic_hist[samplename]["nominal_noweight"].SetDirectory(0)


        ###now plotting systematics
        for sys in self.systematiclist:
            self.sample_systematic_hist[samplename][sys] = {}

            if sys == "QCDscale":
                faileddatasets =  ["WWToLNuQQ_aTGC_13TeV-madgraph-pythia8", "ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1"]
                faileddatasets.append("ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1")
                if (samplename in faileddatasets) or (self.DYdatadriven and samplename.replace("_untagged","") in faileddatasets):
                    self.sample_systematic_hist[samplename][sys]["up"] = None
                    self.sample_systematic_hist[samplename][sys]["down"] = None
                    continue
                plotdir = "QCDScale_intermedianplots/"
                #QCDscalehist = self.GetQCDScaleHist(filename)
                plotname = os.path.join(plotdir, samplename+"_"+channel)
                weight_cuts = "("+ self.cuts +"&&"+ self.channelcuts[channel] +")*"+weight
                hist_up, hist_down = self.QCDScaleSystematic2D(chain, weight_cuts, plotname)
                hist_up.SetName(key+"_"+channel+"_%s"%samplename+"_QCDscaleup")
                hist_down.SetName(key+"_"+channel+"_%s"%samplename+"_QCDscaledown")
                #print "samplename ",samplename, " up ",hist_up.Print("ALL")," down ",hist_down.Print("ALL")," nominal ",self.sample_systematic_hist[samplename]["nominal"].Print("ALL")
                self.sample_systematic_hist[samplename][sys]["up"] = hist_up
                self.sample_systematic_hist[samplename][sys]["down"] = hist_down
                continue

            if sys == "MC_statistical":## question ???
                hist_up = ROOT.TH2F(key+"_"+channel+"_%s"%samplename+"_"+sys+"_up", key+"_"+channel+"_%s"%samplename+"_"+sys+"_up", len(self.xbins)-1, self.xbins,  len(self.ybins) -1, self.ybins)
                hist_down = ROOT.TH2F(key+"_"+channel+"_%s"%samplename+"_"+sys+"_down", key+"_"+channel+"_%s"%samplename+"_"+sys+"_down", len(self.xbins)-1, self.xbins,  len(self.ybins) -1, self.ybins)
                for bin in xrange(self.sample_systematic_hist[samplename]["nominal"].GetNbinsX()):
                    for ybin in xrange(self.sample_systematic_hist[samplename]["nominal"].GetNbinsY()):
                        binvalue = self.sample_systematic_hist[samplename]["nominal"].GetBinContent(bin+1, ybin+1)
                        rawbinvalue = self.sample_systematic_hist[samplename]["nominal_noweight"].GetBinContent(bin+1, ybin+1)
                        weight_thisbin = 1.0
                        if rawbinvalue>0:
                            weight_thisbin = binvalue/rawbinvalue
                        err = sqrt(abs(rawbinvalue))*weight_thisbin
                        hist_up.SetBinContent(bin+1, ybin+1, binvalue + err)
                        hist_down.SetBinContent(bin+1, ybin+1, binvalue - err)
                        if err >  0.1*abs(binvalue) and (key == "TT"):
                            print "warning, in this bin ",bin+1," err ",err, " binvalue ",binvalue, " is larger than 10% of binvalue in TTbar!!!"
                        if bin == 1:
                            print(key+"_"+channel+"_%s"%samplename," binvalue ", binvalue," binvalue_raw ", rawbinvalue," bin_weighthist ", weight_thisbin, " staterr ", err)

                self.sample_systematic_hist[samplename][sys]["up"] = hist_up
                self.sample_systematic_hist[samplename][sys]["down"] = hist_down
                continue

            for plottype in ["up","down"]:

                suffix = sys+"_"+plottype
                hist = ROOT.TH2F(key+"_"+channel+"_%s"%samplename+"_%s"%(suffix), key+"_"+channel+"_%s"%samplename+"_%s"%(suffix), len(self.xbins)-1, self.xbins,  len(self.ybins) -1, self.ybins)
                self.sample_systematic_hist[samplename][sys][plottype] = hist
                
                sysweight = self.SysDict[sys]["weight"].replace("nominal", plottype)
                thisweight = sysweight+"*"+weight
                finalcut = "(" + self.cuts +"&&"+ self.channelcuts[channel] + ")*"+thisweight
                #chlist[samplename].Draw(todraw + ">> " + self.sample_systematic_hist[samplename][sys][plottype].GetName(), finalcut)
                #print "sample ",samplename," file ",self.filedict[key][samplename]['path']," weight ",thisweight," todraw2D ",self.todraw2D," finalcut ",finalcut
                chain.Draw(self.todraw2D + ">> " + self.sample_systematic_hist[samplename][sys][plottype].GetName(), finalcut)
                self.sample_systematic_hist[samplename][sys][plottype].SetDirectory(0)
                

    #################################################################################################
    ##plot systematcs of all sub-processes in one process like TT, sT, DY...  2D
    #################################################################################################
    def plotSystematics2D(self, key, channel):
        
        self.shortname = key
        self.hist_data_untagged = None
        self.sample_systematic_hist = {}
       
        #print "plotSystematics ",key, " channel ", channel
        if self.shortname == "DY" and self.DYdatadriven and (channel == "MuMu" or channel == "ElEl"):
            ##plot data
            self.hist_data_untagged = ROOT.TH2F("data_untagged_"+channel, "data_untagged_"+channel, len(self.xbins)-1, self.xbins,  len(self.ybins) -1, self.ybins)
            Mbtag_weight = "dy_Mbtag_weight"
            untagged_suffix = "_untagged"
            data_untagged_name = "DoubleMuon"
            if channel == "ElEl":
                data_untagged_name = "DoubleEG"
            datafile_untagged =  self.filedict["Data"+untagged_suffix][data_untagged_name]['path']
            ch_d = ROOT.TChain(self.treename)
            ch_d.AddFile(datafile_untagged)
            cut_data = "("+ self.cuts +"&&"+ self.channelcuts[channel] +")*"+Mbtag_weight
            ch_d.Draw(self.todraw2D + ">> " + self.hist_data_untagged.GetName(), cut_data)
            self.hist_data_untagged.SetDirectory(0)

            weight_datadriven = self.nominalWeight+"* "+ Mbtag_weight
            for key_datadriven in self.shortnames_datadriven:
                for samplename_datadriven in self.filedict[key_datadriven].keys():
                    self.plotSystematics_singleFile2D(key, channel, samplename_datadriven, self.filedict[key_datadriven][samplename_datadriven]['path'], weight_datadriven)
        else:
            for iname, samplename in enumerate(self.filedict[self.shortname].keys()):
                self.plotSystematics_singleFile2D(self.shortname, channel, samplename, self.filedict[self.shortname][samplename]['path'], self.nominalWeight)

    
    #################################################################################################
    ##combine subprocess systematcs for one process: like TT, sT, ttV, TT_untagged, sT_untagged, 2D
    #################################################################################################
    def combineSystematics2D(self, shortname, channel):
        ##combine to final plots
        self.shortname = shortname
        suffix = "nominal"
        #print "combineSystematics2D ", shortname, " channel ", channel
        self.finalhist["nominal"] = ROOT.TH2F(self.shortname+"_"+channel, self.shortname+"_"+channel, len(self.xbins)-1, self.xbins,  len(self.ybins) -1, self.ybins)
        allsamples_combine = None
        if "untagged" in self.shortname and self.DYdatadriven and (channel == "MuMu" or channel == "ElEl"):
            allsamples_combine = self.filedict[shortname]
        else:
            allsamples_combine = self.filedict[self.shortname].keys()

        for iname, samplename in enumerate(allsamples_combine):
            self.finalhist["nominal"].Add( self.sample_systematic_hist[samplename]["nominal"] )
        self.finalhist["nominal"].SetDirectory(0)

        rateshift = 1.0
        if self.addRateShift:
            rateshift = self.getTotalRateShift(shortname, channel)

        dummyplots = {}
        totalSys = {}
        for sys in self.systematiclist:
            self.systematic_hist[sys] = {}
            totalSys[sys] = {}
        for plottype in ["up","down"]:

            dummyplots[plottype] = self.cloneDummyHistogram2D( self.finalhist["nominal"] )
            dummyplots[plottype].SetName("dummy"+plottype)
            for sys in self.systematiclist:
                totalSys[sys][plottype] = 0.0
                suffix = "_"+sys+"_"+plottype
                #hist = ROOT.TH2F(self.shortname+"_"+channel+suffix, self.shortname+"_"+channel+suffix, len(self.xbins)-1, self.xbins,  len(self.ybins) -1, self.ybins)
                self.systematic_hist[sys][plottype] = self.cloneDummyHistogram2D( self.finalhist["nominal"] )
                self.systematic_hist[sys][plottype].SetName(self.shortname+"_"+channel+suffix)
                #print "systematic hist name ",self.shortname+"_"+channel+suffix
                for iname, samplename in enumerate(allsamples_combine):
                    if self.sample_systematic_hist[samplename][sys][plottype]:
                        self.systematic_hist[sys][plottype].Add( self.sample_systematic_hist[samplename][sys][plottype] )
                    else:
                        ##if sys shift not available, still need to add the nominal 
                        self.systematic_hist[sys][plottype].Add( self.sample_systematic_hist[samplename]["nominal"] )
                #print "plottype ",plottype, self.systematic_hist[sys][plottype].Print("ALL")," nominal ",self.finalhist["nominal"].Print("ALL")
                self.systematic_hist[sys][plottype].SetDirectory(0)
                for bin in xrange(self.systematic_hist[sys][plottype].GetNbinsX()):
                    for ybin in xrange(self.systematic_hist[sys][plottype].GetNbinsY()):
                        sysvalue = self.systematic_hist[sys][plottype].GetBinContent(bin+1, ybin+1) - self.finalhist["nominal"].GetBinContent(bin+1, ybin+1) ##error: hist_up - hist_nominal
                        #print "plottype ",plottype, " bin ",bin," systematic_hist ",self.systematic_hist[sys][plottype].GetBinContent(bin+1)," nominal ",self.finalhist["nominal"].GetBinContent(bin+1)  ," sysvalue ",sysvalue ," dummyplot bin content ", dummyplots[plottype].GetBinContent(bin+1)
                        totalSys[sys][plottype] = totalSys[sys][plottype] + sysvalue
                        bincontent = dummyplots[plottype].GetBinContent(bin+1, ybin+1) + sysvalue * sysvalue
                        dummyplots[plottype].SetBinContent(bin+1, ybin+1, bincontent)
                totalSys[sys][plottype] = abs(totalSys[sys][plottype])/self.finalhist["nominal"].Integral()

            for bin in xrange( dummyplots[plottype].GetNbinsX()):
                for ybin in xrange( dummyplots[plottype].GetNbinsY()):
                    bincontent = dummyplots[plottype].GetBinContent(bin+1, ybin+1) 
                    nominal =  self.finalhist["nominal"].GetBinContent(bin+1, ybin+1)
                    finalError2 = bincontent + nominal*nominal*(rateshift-1.0)*(rateshift-1.0)
                    dummyplots[plottype].SetBinContent(bin+1, ybin+1, sqrt(finalError2))


            self.finalhist[plottype] = self.finalhist["nominal"].Clone()
            self.finalhist[plottype].SetName(self.shortname+"_"+channel+plottype)
            addweight = 1.0
            #print "plottype ", plottype, " dummpyplot ",dummyplots[plottype].Print("ALL")
            if plottype == "down":
                addweight = -1.0
            self.finalhist[plottype].Add(dummyplots[plottype], addweight)
            self.finalhist[plottype].SetDirectory(0)


        #print "nominal ", self.finalhist["nominal"].Print("ALL")
        ##add one side 
        self.finalhist["one_sided"] = self.cloneDummyHistogram2D( self.finalhist["nominal"] )
        self.finalhist["one_sided"].SetName(self.shortname+"_"+channel+"one_sided")
        self.finalhist["nominal_allSys"] = self.cloneDummyHistogram2D( self.finalhist["nominal"] )
        self.finalhist["nominal_allSys"].SetName(self.shortname+"_"+channel+"allSys")
        
        for bin in xrange(self.finalhist["one_sided"].GetNbinsX()):
            for ybin in xrange(self.finalhist["one_sided"].GetNbinsY()):
                err_up = self.finalhist["up"].GetBinContent(bin + 1, ybin+1)
                err_down = self.finalhist["down"].GetBinContent(bin + 1, ybin+1)
                nominal = self.finalhist["nominal"].GetBinContent(bin + 1, ybin+1)
                staterr_up = self.systematic_hist["MC_statistical"]["up"].GetBinContent(bin+1, ybin+1)
                staterr = abs(staterr_up-nominal) ## should be symetric

                onesided = (abs(nominal - err_up) + abs(nominal - err_down))/2.0
                #print "err_up ",err_up, " err_down ",err_down, " nominal ",nominal," oneside ",onesided," error/nominal ", onesided/nominal
                #sys_stat_total = sqrt(onesided*onesided + nominal)
                self.finalhist["nominal_allSys"].SetBinContent(bin+1, ybin+1, nominal)
                self.finalhist["one_sided"].SetBinContent(bin+1, ybin+1, nominal + onesided)
                self.finalhist["nominal_allSys"].SetBinError(bin+1, ybin+1, staterr)
                self.finalhist["nominal"].SetBinError(bin+1, ybin+1, onesided)


        self.systematic_hist["nominal"] = self.finalhist["nominal"].Clone()
        #print "Integral nominal ", self.finalhist["nominal"].Integral(), " up ",self.finalhist["up"].Integral()," down ",self.finalhist["down"].Integral(), " one sided ",self.finalhist["one_sided"].Integral()
        #print "totalSys ",totalSys," self.systematic_hist ",self.systematic_hist
        for sys in self.systematiclist:
            totalSys[sys]["oneside"] = (totalSys[sys]["up"] + totalSys[sys]["down"])/2.0
            print "Systematic ",sys," samplename ",self.shortname, " up ",totalSys[sys]["up"], " down " ,totalSys[sys]["down"]," oneside ",totalSys[sys]["oneside"]
        print "Final ", self.finalhist["nominal"].Integral(), " one_sided ", self.finalhist["one_sided"].Integral()," error ", self.finalhist["one_sided"].Integral()/self.finalhist["nominal"].Integral()-1.0


    #################################################################################################
    ###  summerize systematics 2D
    #################################################################################################
    def runSystematics2D(self, shortname, channel):
        self.plotSystematics2D(shortname, channel)
        if shortname == "DY" and self.DYdatadriven and (channel == "MuMu" or channel == "ElEl"):
            datadriven_finalhist ={}
            for name in self.shortnames_datadriven:
                self.combineSystematics2D(name, channel)
                datadriven_finalhist[name] = self.finalhist.copy()
                self.channel_shortname_systematic_hist[channel][name] = self.systematic_hist.copy()
                #print "Data-driven ", datadriven_finalhist[name]["nominal"].Print("ALL")

            ### calculate final hist for data-driven DY    
            self.finalhist = {}
            for plottype in ["nominal","up","down","one_sided"]:
                suffix = "_"+plottype
                if plottype == "nominal":
                    suffix = ""
                self.finalhist[plottype] = self.cloneDummyHistogram2D(self.hist_data_untagged)
                self.finalhist[plottype].SetName(shortname + "_"+ channel+ suffix)

            for  bin in xrange(self.hist_data_untagged.GetNbinsX()):
                for  ybin in xrange(self.hist_data_untagged.GetNbinsY()):
                    err_up_square = 0.0
                    err_down_square = 0.0
                    nominal_mc = 0.0
                    data_binvalue = self.hist_data_untagged.GetBinContent(bin+1, ybin+1)
                    for name in self.shortnames_datadriven:
                        thisnominal = datadriven_finalhist[name]["nominal"].GetBinContent(bin+1, ybin+1)
                        thisup = datadriven_finalhist[name]["up"].GetBinContent(bin+1, ybin+1)
                        thisdown = datadriven_finalhist[name]["down"].GetBinContent(bin+1, ybin+1)
                        err_up_square = err_up_square + (thisnominal-thisup)*(thisnominal-thisup)
                        err_down_square = err_down_square + (thisnominal-thisdown)*(thisnominal-thisdown)
                        nominal_mc  = nominal_mc+thisnominal
                    one_sided_err = (sqrt(err_up_square)+sqrt(err_down_square))/2.0
                    self.finalhist["nominal"].SetBinContent(bin+1, ybin+1, data_binvalue - nominal_mc)
                    self.finalhist["up"].SetBinContent(bin+1,ybin+1,  data_binvalue - nominal_mc + sqrt(err_up_square))
                    self.finalhist["down"].SetBinContent(bin+1, ybin+1, data_binvalue - nominal_mc - sqrt(err_down_square))
                    self.finalhist["one_sided"].SetBinContent(bin+1, ybin+1, data_binvalue - nominal_mc + one_sided_err)
                    self.finalhist["nominal"].SetBinError(bin+1, ybin+1, one_sided_err)
            #print "Data-Driven DY final ", self.finalhist["nominal"].Print("ALL")


        else:
            self.combineSystematics2D(shortname, channel)
            self.channel_shortname_systematic_hist[channel][shortname] = self.systematic_hist.copy()
        #print "shortname ",shortname, " all process systematics ", self.channel_shortname_systematic_hist

