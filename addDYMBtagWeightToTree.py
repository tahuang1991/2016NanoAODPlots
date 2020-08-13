import ROOT
import argparse
import os
import numpy as np
from array import array


parser = argparse.ArgumentParser(description='Compute event category fraction on a given sample')
parser.add_argument('-i', '--inputfile',dest= "inputfile", type=str, metavar='STR', help='input file')
options = parser.parse_args()


flavours = ['b','c','l']
def getEff1D(teff, x):
    hist = teff.GetCopyPassedHisto()
    xbin = hist.GetXaxis().FindBin(x)
    return teff.GetEfficiency(xbin)

def getEff2D(teff, x, y):
    hist = teff.GetCopyPassedHisto()
    xbin = hist.GetXaxis().FindBin(x)
    ybin = hist.GetYaxis().FindBin(y)
    bin = hist.GetBin(xbin, ybin)
    return teff.GetEfficiency(bin)

def get_Mbtag_weight(bdt_value,  pt1, eta1,  pt2, eta2, allfracs, Teffs):
    weight = 0.0
    for f1 in flavours:
        for f2 in flavours:
            thisfrac = allfracs[f1+f2]
            eff1 = getEff2D(Teffs[f1], pt1, eta1)
            eff2 = getEff2D(Teffs[f2], pt2, eta2)
	    frac = getEff1D(thisfrac, bdt_value)
	    weight = weight + frac*eff1*eff2
            #print "flavours ",f1," ",f2," eff1 ",eff1, " eff2 ",eff2 ," frac ", frac," current weight ", weight, " frac*eff1*eff2 ",frac*eff1*eff2
    return weight 

dy_frac_eff_f = "dy_frac_eff_combined.root"
dytf = ROOT.TFile(dy_frac_eff_f, "READ")
allfractions = {}
allTeffs = {}

for f1 in flavours:
    for f2 in flavours:
        allfractions[f1+f2] = dytf.Get(f1+f2+"_frac")
        allfractions[f1+f2].SetDirectory(0)
for f in flavours:
    eff = None
    if f != 'l':
       eff = dytf.Get("btagging_eff_on_"+f)
    else:
       eff = dytf.Get("mistagging_eff_on_light")
    allTeffs[f] = eff
    allTeffs[f].SetDirectory(0)
         
dytf.Close()
#print "allfractions ",allfractions," alleff ", allTeffs



print "=============================================================="
print "inputfile ", options.inputfile
print "=============================================================="

_rootBranchType2PythonArray = { 'b':'B', 'B':'b', 'i':'I', 'I':'i', 'F':'f', 'D':'d', 'l':'L', 'L':'l', 'O':'B' }

treename = "Friends"
tf = ROOT.TFile(options.inputfile, "update")
chain = tf.Get( treename )


rootBranchType = "F"
bdt_value =  array(_rootBranchType2PythonArray[rootBranchType], [0.0])
#dy_Mbtag_weight =  array(_rootBranchType2PythonArray[rootBranchType], [0.0])
dy_Mbtag_weight_v2 =  array(_rootBranchType2PythonArray[rootBranchType], [0.0])
#br_bdt_value = chain.Branch("bdt_value", bdt_value ,"bdt_value/%s" % (rootBranchType))
#br_dy_Mbtag_weight = chain.Branch("dy_Mbtag_weight", dy_Mbtag_weight , "dy_Mbtag_weight/%s" % (rootBranchType))
br_dy_Mbtag_weight_v2 = chain.Branch("dy_Mbtag_weight_v2", dy_Mbtag_weight_v2, "dy_Mbtag_weight_v2/%s" % (rootBranchType))
def initBranchVals():
    #bdt_value[0] = -99
    #dy_Mbtag_weight[0] = 0.0
    dy_Mbtag_weight_v2[0] = 0.0

#def fillBranches():
#    print "fillBranches():  bdt_value ",bdt_value[0]," weight ",dy_Mbtag_weight[0]
#    br_bdt_value.Fill()
#    br_dy_Mbtag_weight.Fill()

entries = None
if not entries:
    entries = chain.GetEntries()



chain.SetBranchStatus("*", 0)

chain.SetBranchStatus("is*", 1)
#chain.SetBranchStatus("lep*", 1)
chain.SetBranchStatus("jet*", 1)
chain.SetBranchStatus("ll*", 1)
chain.SetBranchStatus("jj*", 1)
chain.SetBranchStatus("jj_pt", 1)
chain.SetBranchStatus("ht*", 1)
chain.SetBranchStatus("nJetsL*", 1)
chain.SetBranchStatus("*DPhi*", 1)
#chain.SetBranchStatus("*partonFlavour*", 1)
chain.SetBranchStatus("*weight*", 1)
chain.SetBranchStatus("cosThetaStar", 1)
chain.SetBranchStatus("sample_weight", 1)
chain.SetBranchStatus("bdt_value", 1)
chain.SetBranchStatus("dy_Mbtag_weight_v2", 1)

"""
bdt_tmva_variables = [
        "jet1_pt",
        "jet1_eta",
        "jet2_pt",
        "jet2_eta",
        "jj_pt",
        "ll_pt",
        "ll_eta",
        "llmetjj_DPhi_ll_met",
        "ht",
        "nJetsL"

]
#bdt_label = "2016_12_18_BDTDY_bb_cc_vs_rest_7var_ht_nJets"
#bdt_label = "2017_02_17_BDTDY_bb_cc_vs_rest_10var"
bdt_label = "2018_04_12_BDTDY_bb_cc_vs_rest_10var"
#FIXME
bdt_xml_file = "DYBDTTraining/weights/{}_kBDT.weights.xml".format(bdt_label)
print "bdt_xml_file ",bdt_xml_file

dict_tmva_variables = { var: array('f', [0]) for var in bdt_tmva_variables }
m_reader = ROOT.TMVA.Reader("Silent=1")
for var in bdt_tmva_variables:
    m_reader.AddVariable(var, dict_tmva_variables[var])
m_reader.BookMVA(bdt_label, bdt_xml_file)
entries = None
print("Loading chain...")
if not entries:
    entries = chain.GetEntries()
print("Done.")
"""
print("Adding medium btagging weight for  %d events." % entries)


for i in range(0, entries):
    chain.GetEntry(i)

    initBranchVals()
    if (i % 10000 == 0):
        print("Event %d over %d" % (i + 1, entries))


    if not (chain.isElEl or chain.isMuMu or chain.ll_M>12):
        #print "ignore Muel channel"
        #fillBranches()
        continue
        br_bdt_value.Fill()
        br_dy_Mbtag_weight.Fill()


    def get_value(object, val):
	return getattr(object, val)

    #for var in bdt_tmva_variables:
    #    # Special treatment for variables not retrieved from the base object
    #    dict_tmva_variables[var][0] = get_value(chain, var)

    #bdt_value[0] = m_reader.EvaluateMVA(bdt_label)
    bdt_value[0] = chain.bdt_value
    dy_Mbtag_weight_v2[0] = get_Mbtag_weight(bdt_value[0], chain.jet1_pt, abs(chain.jet1_eta), chain.jet2_pt, abs(chain.jet2_eta), allfractions, allTeffs)
    #print " bdt_value ",bdt_value[0]," weight ",dy_Mbtag_weight_v2[0]
#    #fillBranches()
#    br_bdt_value.Fill()
    br_dy_Mbtag_weight_v2.Fill()
#    
#
#   
#
chain.Write()
tf.Close()
