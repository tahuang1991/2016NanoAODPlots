import numpy as np
import ROOT
import sys
sys.argv.append('-b')
import matplotlib.pyplot as plt
import os

def f(t):
    return np.exp(-t) * np.cos(2*np.pi*t)

masspoints = np.asarray([260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900])
MTonly = {800: {'isMuMu': 70.95057123028768, 'isElEl': 61.577868874018584, 'isElMu': 67.53708922016595}, 450: {'isMuMu': 48.49919040620318, 'isElEl': 39.1116895624184, 'isElMu': 45.85987819858824}, 260: {'isMuMu': 29.545854031119784, 'isElEl': 19.24653787287462, 'isElMu':     26.039224591551832}, 550: {'isMuMu': 60.78347353280954, 'isElEl': 49.484077132048206, 'isElMu': 55.362923659545864}, 650: {'isMuMu': 67.52311194349906 , 'isElEl': 56.25571235467507, 'isElMu': 61.95240667781626}, 300: {'isMuMu': 28.45511365986146, 'isElEl': 19.801793984157683, 'isElMu': 25.147529482518383}, 270: {'isMuMu': 28.246059454189822, 'isElEl': 18.617255473983548, 'isElMu': 25.43536357474379}, 400: {'isMuMu': 41.336034387388814, 'isElEl': 31.729708505350754, 'isElMu': 37.84382853317458}, 500: {'isMuMu': 55.03049202688094, 'isElEl': 45.22074949053289, 'isElMu': 51.518357060406245}, 750: {    'isMuMu': 68.91683051592886, 'isElEl': 59.52602477699557, 'isElMu': 65.78754189387902}, 600: {'isMuMu': 62.941651213847365, 'isElEl': 53.026766869769496, 'isElMu': 60.16692021097901}, 900: {'isMuMu': 72.10203660139918, 'isElEl': 62.974648267519484, 'isElMu': 69.01881301454304}, 350: {'isMuMu': 32.99032066963954, 'isElEl': 24.538713030935376, 'isElMu': 29.549837199581567}}
MTandMT2 = {800: {'isMuMu': 71.68697360146099, 'isElEl': 61.94303692562975, 'isElMu': 68.16428699876205}, 450: {'isMuMu': 51.30493359614676, 'isElEl': 41.18494126435214, 'isElMu': 48.18695569492916}, 260: {'isMuMu': 31.55648829862514, 'isElEl': 20.75878196090403, 'isElMu': 28.025593132879536}, 550: {'isMuMu': 63.05993703634365, 'isElEl': 51.398253150811975, 'isElMu': 57.76609709738885}, 650: {'isMuMu': 69.20277461396681, 'isElEl': 57.331083371362226, 'isElMu': 63.37946786580976}, 300: {'isMuMu': 29.19658267628301, 'isElEl': 20.53072096774969, 'isElMu': 25.951511195476773}, 270: {'isMuMu': 29.860215216164818, 'isElEl': 19.81367506176357, 'isElMu': 27.264840744147143}, 400: {'isMuMu': 43.338683080725126, 'isElEl': 33.345839488670954, 'isElMu': 39.469420326703734}, 500: {'isMuMu': 57.71442164151281, 'isElEl': 47.43823330535144, 'isElMu': 54.04243075752646}, 750: {'isMuMu': 69.96836945292924, 'isElEl': 60.198651923859636, 'isElMu': 66.63548848470904}, 600: {'isMuMu': 65.05255427853517, 'isElEl': 54.384822944936154, 'isElMu': 61.89273844725654}, 900: {'isMuMu': 72.49390424304302, 'isElEl': 63.11034731301055, 'isElMu': 69.3143869295379}, 350 : {'isMuMu': 33.64776295461901, 'isElEl': 24.945041545432858, 'isElMu': 30.047378152011106}}
MTandMT2_MJJ =  {800: {'isMuMu': 71.78775852223859, 'isElEl': 62.05355805274861, 'isElMu': 68.30755942586815}, 450: {'isMuMu': 52.20520847761798, 'isElEl': 41.85226020047538, 'isElMu': 49.137507725814714}, 260: {'isMuMu': 32.05905852595361, 'isElEl': 20.959103330335623, 'isElMu': 28.206537259318246}, 550: {'isMuMu': 63.53069481184263, 'isElEl': 51.59604224748493, 'isElMu': 58.12911111813868}, 650: {'isMuMu': 69.31845536148275, 'isElEl': 57.509197437372514, 'isElMu': 63.543916311615334}, 300: {'isMuMu': 30.052740380604266, 'isElEl': 20.87692485172653, 'isElMu': 26.539694307287256}, 270 : {'isMuMu': 30.394630973036726, 'isElEl': 20.28571993611892, 'isElMu': 27.455655108008298}, 400: {'isMuMu': 44.74482847476057, 'isElEl': 34.43361585324694, 'isElMu': 40.70317693933215}, 500: {'isMuMu': 58.32251649252591, 'isElEl': 47.77172779763888, 'isElMu': 54.54764148434456}, 750: {'isMuMu': 70.10893926678736, 'isElEl': 60.24952949389543, 'isElMu': 66.72351111274993}, 600: {'isMuMu': 65.2525173383583, 'isElEl': 54.565939154125424, 'isElMu': 62.10985187742745}, 900: {'isMuMu': 72.65824592903897, 'isElEl': 63.297918193600836, 'isElMu': 69.44872607802813}, 350: {'isMuMu': 35.43766952089536, 'isElEl': 25.97656404526518, 'isElMu': 31.169519726343612}}
#sb1 = np.asarray([48.74, 46.96, 47.62, 56.77, 72.51, 87.553, 99.24, 108.93, 115.81, 121.91, 127.60, 131.59, 134.45])
#sb2 = np.asarray([52.41, 49.96, 49.02, 57.8, 75.88, 92.26, 104.19, 113.32, 119.25, 124.67, 129.26, 132.66, 134.98])
#sb3 = np.asarray([53.04, 50.70, 50.31, 60.42,78.34, 93.87, 105.14, 114.01, 119.61, 124.93, 129.48, 132.91, 135.28])

def getSoverB_mass(MTonly, masslist, channel):
    sb = []
    for mass in masslist:
        sb.append(MTonly[int(mass)][channel])
    #print "channel ", channel, " s over b ",sb, len(sb), " len of masspoints ",len(masspoints)
    return sb



outdir = "MT_MT2_MJJ_SoverB_20180214_bin80/"
os.system("mkdir -p "+outdir)

colors = ['blue', 'green', 'red','orange', 'black']
for channel in ['isMuMu', 'isElMu','isElEl']:
    sb1 =  np.asarray( getSoverB_mass(MTonly, masspoints, channel) )
    sb2 =  np.asarray( getSoverB_mass(MTandMT2, masspoints, channel) )
    sb3 =  np.asarray( getSoverB_mass(MTandMT2_MJJ, masspoints, channel) )
    ratio1 = np.divide(sb2, sb1)
    ratio2 = np.divide(sb3, sb1)
    #print "sb1 ",sb1," masspoints ",masspoints 
#fig  = plt.figure(1)
#fig, (ax1,ax2) = plt.subplots(1, 2)
    #fig, ax1 = plt.subplots(1, 1)
    fig, (ax1, ax2) = plt.subplots(2,1, gridspec_kw = {'height_ratios':[2, 1]})
    ax1.plot(masspoints, sb1, color=colors[0], lw=3, label= "MTonly")
    
    ax1.plot(masspoints, sb2, color=colors[1], lw=3, label= "MT+MT2")
    ax1.plot(masspoints, sb3, color=colors[2], lw=3, label= "MT+MT2+Mjj")
#label.set_color('k')
    ax1.margins(0.05)
#ax1.set_xlim(xlim);
#ax1.set_ylim(ylim);
    #ax1.set_xlabel("Mass")
    ax1.set_ylabel("$S/\sqrt{B}$")
#f1ig.set_tight_layout(True)
    str_M = "Sensitivity estimation(No systematics) for all mass points, %s channel"%channel[2:]
    ax1.grid(color='k', linestyle='--', linewidth=1)
    axlegs = ax1.legend(loc='center right', numpoints=1, title =" ", frameon=False, prop={'size': 12})
#axlegs.get_title().set_fontsize('14')
#axlegs.set_title(str_M, {'size' : 14})
    ax1.set_title(str_M)
    for color,text in zip(colors,axlegs.get_texts()):
        text.set_color(color)

    ax2.plot(masspoints, ratio1, color=colors[1], lw=3, linestyle='dashed', label= "MT+MT2")
    ax2.plot(masspoints, ratio2, color=colors[2], lw=3, linestyle='dashed', label= "MT+MT2+Mjj")
    ax2legs = ax2.legend(loc='center right', numpoints=1, title =" ", frameon=False, prop={'size': 12})
    for color,text in zip(colors[1:],ax2legs.get_texts()):
        text.set_color(color)
    ax2.grid(color='k', linestyle='--', linewidth=1)
    ax2.set_xlabel("Mass [GeV]")
    ax2.set_ylabel("Ratio: over MTonly")
#plt.subplot(212)
#ax2.plot(t2, np.cos(2*np.pi*t2), 'r--', label="ax2label")
#leg = ax2.legend(loc='center right', numpoints=1, frameon=False)
#h, l = ax2.get_legend_handles_labels()
#leg.texts[0].set_color('r')
#print " leg.texts ",leg.texts
#print " ax2 label ",ax2.get_label().set_color("red")
#plt.show()
#fig.savefig("twoplotinonefigure.png")
    plt.subplots_adjust( hspace=0.001)
    #fig.savefig(outdir+"SoverBsummary_3models_%s_test.pdf"%channel)

wsfile = ROOT.TFile("","READ")
