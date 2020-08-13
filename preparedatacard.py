import os
import ROOT
import datetime

def generatedatacard(outfile, rootfile, channels, processes, mass, addobservation):
    datacardfile = open(outfile, "write")
    imax = len(channels)
    jmax = len(processes)
    datacardfile.write("imax %d\n"%imax)
    datacardfile.write("jmax %d\n"%jmax)
    datacardfile.write("kmax *\n")  ## uncertainty sources
    datacardfile.write("#-----------------------------------------------------------------------------\n")## comment out
    rootfilename = rootfile.split('/')[-1]
    datacardfile.write("shapes * * %s $PROCESS_$CHANNEL_M%d $PROCESS_$SYSTEMATIC\n"%(rootfilename, mass))
    datacardfile.write("#-----------------------------------------------------------------------------\n")## comment out
    ## channels
    datacardfile.write("bin\t\t")
    for chan in channels:
	datacardfile.write(chan+"\t")
    datacardfile.write("\n")
    rfile = ROOT.TFile(rootfile, "READ") ## read rate from histogram
    if addobservation:
        datacardfile.write("observation\t\t")
	for chan in channels:
	    hist = rfile.Get("data_obs_{channel}_M{mass}".format(channel = chan, mass = mass))
	    datacardfile.write("{rate}\t".format(rate = hist.Integral()))
	datacardfile.write("\n")
    else:
        datacardfile.write("#observation, not available yet\n")
    datacardfile.write("#-----------------------------------------------------------------------------\n")## comment out
    ##start to set rate
    nratelines = 4
    ratelineheads = {"bin","process","process","rate"}
    for iline,head in enumerate(ratelineheads):
        datacardfile.write("{head}\t\t".format(head = head))
        for i, chan in enumerate(channels):
	    for j, process in enumerate(processes): 
                if iline == 0:
                    datacardfile.write("{channel}\t".format(channel = chan))
                elif iline == 1:
                    datacardfile.write("{process}\t".format(process = process))
                elif iline == 2:
                    datacardfile.write("%d\t"%(j))
                elif iline == 3:
                    hist = rfile.Get("{process}_{channel}_M{mass}".format(process = process, channel = chan, mass = mass))
                    datacardfile.write("%.1f"%(hist.Integral()))
        datacardfile.write("\n")
    datacardfile.write("#-----------------------------------------------------------------------------\n")## comment out
    ## add uncertainty, should be done in a smarter way?
    funcertainty = open("allUncertainties_1.txt","read")
    for line in funcertainty:
        datacardfile.write(line)
     
    rfile.Close()
    datacardfile.close()

def alldatacards(rootfile, masspoints, addobservation, outdir):
    os.system("cp %s %s"%(rootfile, outdir))
    channels = ["MuMu","MuEl","ElEl"]
    processes = ["signal","TT","DY","sT"]
    for mass in masspoints:
        outfile = os.path.join(outdir, "shape_M%d.txt"%mass)
        generatedatacard(outfile, rootfile, channels, processes, mass, addobservation)

suffix = 'addNN'
output_suffix = '{:%Y-%m-%d}_{}'.format(datetime.date.today(), suffix)
datacarddir = "Datacards/"
os.system("mkdir -p "+datacarddir)
outdir = os.path.join(datacarddir, output_suffix)
os.system("mkdir -p "+outdir)
rootfile = "HistForLimits/Hhh_FinalBGYield_5pbsignalindata.root"
masslist = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900]
alldatacards(rootfile, masslist, False, outdir)
