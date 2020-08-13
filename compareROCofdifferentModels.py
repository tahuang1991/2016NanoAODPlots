import os
import plotTools
from common import *
import keras

from sklearn import metrics
import matplotlib.pyplot as plt

import ROOT
inputs = [
        "jj_pt", 
        "ll_pt",
        "ll_M",
        "ll_DR_l_l",
        "jj_DR_j_j",
        "llmetjj_DPhi_ll_jj",
        "llmetjj_minDR_l_j",
        "llmetjj_MTformula",
        "isSF"
        ]
HMEinputs = [
	"hme_h2mass_reco",
	#"hme_mean_weight2_reco",
	"hme_stddev_reco",
	"hme_entries_reco",
	"hme_entry_peak_reco",
	"hme_mostprob_offshellWmass_reco",
	"hme_stddev_offshellWmass_reco"

	   ]
HMEinputs2 = [
	"hme_h2mass_reco",
	"hme_stddev_reco",
	"hme_entries_reco",
	"hme_h2mass_weight2_reco",
	"hme_stddev_weight2_reco",
	"hme_entries_weight2_reco"

	   ]

inputsnoHME = inputs
inputswithHME = inputs+HMEinputs
inputswithHME2 = inputs+HMEinputs2


cut = "(91 - ll_M) > 15"
cut2 = cut +"&& hme_h2mass_reco>=250.0"

allinputs = [inputsnoHME, inputswithHME, inputswithHME2]
allcuts = [cut, cut2, cut2]
print "allinputs ",allinputs, " allcuts ",allcuts

models = {
	'noHME': {
				#'file': 'hh_resonant_trained_models/2017-01-24_260_300_400_550_650_800_900_dy_estimation_from_BDT_new_prod_on_GPU_deeper_lr_scheduler_100epochs/hh_resonant_trained_model.h5',
	     'file': '/Users/taohuang/Documents/DiHiggs/20170530/20170905_Louvain_10k/HH_Keras/HHTools/mvaTraining/hh_resonant_trained_models_NewNtuples_HME_v12_noHME/2017-09-20_90epochs/hh_resonant_trained_model.h5',
	     'legend': 'Training without HME',
	     'color': '#1f77b4'
	},
	'withHME2':{
	    'file': '/Users/taohuang/Documents/DiHiggs/20170530/20170905_Louvain_10k/HH_Keras/HHTools/mvaTraining/hh_resonant_trained_models_NewNtuples_HME_v11_weight/2017-09-19_90epochs/hh_resonant_trained_model.h5',
            'legend': 'Training with HME(no weight)',
	    'color': '#1f77b4'
	},
	'withHME1':{
	    'file' : '/Users/taohuang/Documents/DiHiggs/20170530/20170905_Louvain_10k/HH_Keras/HHTools/mvaTraining/hh_resonant_trained_models_NewNtuples_HME_v14_weight_weight2/2017-09-21_90epochs/hh_resonant_trained_model.h5',
            'legend': 'Training with HME(weight from offshell Wmass)',
            'color': '#1f77b4'

		}

        }

print "model keys ",models.keys()," 0 ",models.keys()[0]
#masses = sorted(models['dedicated'].keys())
signal_masses = resonant_signal_masses

#
#dataset = DatasetManager(inputs, resonant_weights, cut)
#dataset.load_resonant_signal(signal_masses, add_mass_column=True, fraction=1)
#dataset.load_backgrounds(add_mass_column=True)
#
#all_model = keras.models.load_model(models['all']['file'])

output_dir = 'CompareROCplots_v3'
output_name = 'roc_comparison_resonant.pdf'

outputrootfile = output_dir +"/" + output_name +".root"
#outputrootfile = output_name +".root"
allXs = []; allYs = []
allaucs = []

fig = plt.figure(1, figsize=(7, 7))
ax = fig.add_subplot(111)


line_styles = ['-', '--', '-.']

index = 0
#for case in models:
for case in ['noHME']:
    
    dataset = DatasetManager(allinputs[index], resonant_weights, allcuts[index])
    dataset.load_resonant_signal(signal_masses, add_mass_column=True, fraction=1)
    dataset.load_backgrounds(add_mass_column=True)
    
    #dataset.load_resonant_signal()
    train_model =  keras.models.load_model(models[case]['file'])
    # First, the super model
    os.system("mkdir -p "+output_dir+"_"+case)
    #draw_resonant_training_plots(train_model, dataset, output_dir+"_"+case, split_by_mass=False):
    all_signal_predictions = dataset.get_signal_predictions(train_model)
    all_background_predictions = dataset.get_background_predictions(train_model)

    all_n_signal, _, binning = plotTools.binDataset(all_signal_predictions, dataset.get_signal_weights(), bins=50, range=[0, 1])
    all_n_background, _, _ = plotTools.binDataset(all_background_predictions, dataset.get_background_weights(), bins=binning)
    x, y = plotTools.get_roc(all_n_signal, all_n_background)
    ax.plot(x, y, line_styles[index], color=models[case]['color'], lw=2, label=models[case]['legend'])
    allXs.append(x)
    allYs.append(y)
    
    auc = metrics.auc(x, y, reorder=True)
    allaucs.append(auc)
    print " model ",models[case]['file']," train_model ",train_model," AUC ",auc," Xs ",x," Y ", y
    

    #print "roc curve x ",x," type (x) ",type(x)
    index += 1

"""
ax.set_xlabel("Background efficiency", fontsize='large')
ax.set_ylabel("Signal efficiency", fontsize='large')

ax.margins(0.05)
fig.set_tight_layout(True)

ax.legend(loc='lower right', numpoints=1, frameon=False)
fig.savefig(os.path.join(output_dir, output_name))

print("Comparison plot saved as %r" % os.path.join(output_dir, output_name))

rfile = ROOT.TFile(outputrootfile,"RECREATE")

c1 = ROOT.TCanvas()
mg = ROOT.TMultiGraph()

colors = [ROOT.kBlack, ROOT.kRed, ROOT.kGreen+2]
leg = ROOT.TLegend(.6, 0.2, 85, .4)
leg.SetFillColor(ROOT.kWhite)
glist = []
for i in range(0, index):
    case = models.keys()[i]
    tp_auc = ROOT.TParameter(float)("auc"+case, auc)
    g = ROOT.TGraph(len(allXs[i]), allXs[i], allYs[i])
    g.SetName(case+"_ROC")
    g.SetLineColor(colors[i])
    g.SetLineWidth(2)
    g.GetXaxis().SetTitle("Background efficiency")
    g.GetYaxis().SetTitle("Signal efficiency")
    leg.AddEntry(g,case,"pl")
    auc = allaucs[i]
    g.SetTitle("ROC curve, Signal efficiency Vs background efficiency, AUC=%.4f"%auc)
    g.Write()
    tp_auc.Write()
    glist.append(g)
    #rfile.Close()
    #g.Print("ALL")
    #mg.Add(g)


for g in glist:
	mg.Add(g)
mg.Draw("AP")
leg.Draw("same")
mg.GetXaxis().SetTitle("Background efficiency")
mg.GetYaxis().SetTitle("Signal efficiency")

mg.SetName("ComparisonHME")
mg.SetTitle("Comparison between including HME and NOT including HME in DNN training")
#mg.GetXaxis().SetTitle("Background efficiency")
#mg.GetYaxis().SetTitle("Signal efficiency")

mg.Write()
c1.SaveAs(output_dir+"/"+"roc_comparison_resonant_root.C")
c1.SaveAs(output_dir+"/"+"roc_comparison_resonant_root.pdf")
rfile.Close()
"""
