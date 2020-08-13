import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cbook as cbook
from sklearn import metrics
import datetime
from common import *
from SampleTypeLUT import *


colors = ['blue', 'green', 'red','orange', 'black']
def overImposePlots(cvspaths, legs, xlabel, ylabel, xlim, ylim, title, output_dir, output_name):
  
  def extractmass(cvsname):
      name = cvsname.split("/")[-1]
      if "fixed_M" in name:
          M = 0
          for x in name.split("_"):
              try:
                  M = float(x)
                  return "Signal: M=%d GeV"%M
              except ValueError:
                  pass
      else:
          return "Signal: all mass points"

#fig, (ax, ax2) = plt.subplots(1, 2)
  filename = os.path.join(output_dir, output_name)
  fig = plt.figure(1, figsize=(7, 7), dpi=300)
  fig2 = plt.figure(2, figsize=(7, 7), dpi=300)
  fig.clear()
  fig2.clear()

  ax = fig.add_subplot(111) # Create an axes instance

  ax2 = fig2.add_subplot(111) # Create an axes instance
  
  str_M = ""
  for i , thiscvs in enumerate(cvspaths):
      if not os.path.isfile(thiscvs[0]):
          return
      x1 =  np.genfromtxt(thiscvs[0], delimiter=',')#background eff
      y1 =  np.genfromtxt(thiscvs[1], delimiter=',')#Signal eff
      str_M = extractmass(thiscvs[0])
      ax.plot(x1, y1, '-', color=colors[i], lw=3, label=legs[i])
      binwidth = 1.0/len(x1)
      bincentres =np.asarray( [ 1.0 - binwidth*(j+0.5) for j in range(len(x1))] )
      #print " bincentres ", bincentres, " len ", len(bincentres) ," type ",type(bincentres) ," x1 ",x1," type(x1) ",type(x1)
      ax2.plot(bincentres, x1, '--', color=colors[i], lw=3)
      ax2.plot(bincentres, y1, '-', color=colors[i], lw=3,  label=legs[i])
      #ax2.get_label().set_color( colors[i] )
      
      
      
  ax.margins(0.05)
  #if iR==0:
  ax.set_xlim(xlim);
  ax.set_ylim(ylim);
  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)
  #fig.set_tight_layout(True)
  ax.grid(color='k', linestyle='--', linewidth=1)
  axlegs = ax.legend(loc='center right', numpoints=1, title =str_M, frameon=False, prop={'size': 12})
  #axlegs.get_title().set_fontsize('14')
  axlegs.set_title(str_M, {'size' : 14})
  for color,text in zip(colors,axlegs.get_texts()):
    text.set_color(color)

  if title != '':
	  ax.set_title(title)
	  ax2.set_title("Signal and Background efficiency, "+title.split(',')[-1])

  #fig.savefig(os.path.join(output_dir, output_name))
  fig.savefig(filename+'.pdf', bbox_inches='tight')

  ax2.margins(0.05)
  ax2.set_xlabel("NN ouput")
  ax2.set_ylabel("Efficiency")
  ax2.grid(color='k', linestyle='--', linewidth=1)
  ax2legs = ax2.legend(loc='center right', numpoints=1,  title =str_M,frameon=False, prop={'size': 12})
  ax2legs.get_title().set_fontsize('14')
  ax2.text(.12,.72 ,"signal: solid,  background: dash", {'color': 'k', 'fontsize': 12})
  for color,text in zip(colors,ax2legs.get_texts()):
    text.set_color(color)

  fig2.set_tight_layout(True)
  #fig2.savefig(os.path.join(output_dir, output_name))
  fig2.savefig(filename+'_EfficiencyVSNNoutput.pdf', bbox_inches='tight')


  plt.close()


def overImposePlots_roc(ModelLUT, modellist, output_folder):
  #overall performance
  xlim = [-0.01, 0.5]
  ylim = [.4, 1.01]
  cvspathlist_overall = []
  legs_overall = []
  for key in modellist:
    cvspath_x = os.path.join(ModelLUT[key]['workingdir'],'roc_curve_X.cvs')  
    cvspath_y = os.path.join(ModelLUT[key]['workingdir'],'roc_curve_Y.cvs')  
    cvspathlist_overall.append([cvspath_x, cvspath_y])
    legs_overall.append(ModelLUT[key]['legend'])

  output_name = 'roc_curve_fixed_M_allmass'
  title = 'ROC curve from DNN, Signal: all mass points'
  overImposePlots(cvspathlist_overall, legs_overall, 'Background efficiency', 'Signal efficiency', xlim, ylim, title, output_folder, output_name)
  for mass in resonant_signal_masses:
  #for mass in [260]:
    cvspathlist = []
    legs = []
    for key in modellist:
         filepath = os.path.join(ModelLUT[key]['workingdir'],'splitted_by_mass') 
         cvspath_x = os.path.join(filepath, 'roc_curve_fixed_M_%d_X.cvs'%mass) 
         cvspath_y = os.path.join(filepath, 'roc_curve_fixed_M_%d_Y.cvs'%mass) 
         cvspathlist.append([cvspath_x, cvspath_y])
         legs.append(ModelLUT[key]['legend'])

    title = 'ROC curve from DNN, Signal: M=%d GeV'%mass
    output_name = 'roc_curve_fixed_M_%d_MT2AndMTcomparison'%mass
    overImposePlots(cvspathlist, legs, 'Background efficiency', 'Signal efficiency', xlim, ylim, title, output_folder, output_name)




suffix = '2018_MT_MT2_MJJ_HME_addgrid'
#suffix = '2018_Test'
#suffix = '2017and2018_MTandMT2_HMEMJJ'
output_suffix = '{:%Y-%m-%d}_{}'.format(datetime.date.today(), suffix)
output_folder = os.path.join('ModelComparison', output_suffix)
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
#overImposePlots_roc(ModelLUT, ['MTandMT2','MTandMT2_MJJ','MTandMT2_HME','MTandMT2_HMEMJJ'], output_folder)
#overImposePlots_roc(ModelLUT, ['MTandMT2_HMEMJJ','MTandMT2_HMEMJJ_new'], output_folder)
overImposePlots_roc(ModelLUT, ['MTonly','MTandMT2','MTandMT2_HME'], output_folder)
