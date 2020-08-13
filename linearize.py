import ROOT

masslist =      [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900]
process = 	["ElEl", "MuEl", "MuMu"]
varlist = 	["MTonly", "MTandMT2", "MTandMT2_MJJ"]

workdir = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy2D_HMEbinsv4_1p15/"
for var in varlist:
    #workdir = var+"/"
  new = ROOT.TFile.Open(workdir+"linear_{var}.root".format(var = var), "recreate")
  new.Close()

  for mass in masslist:
    f = ROOT.TFile.Open(workdir+"Hhh_FinalBGYield_xsec1pb_NNvsHME_nnout_{var}_nnstep0p04_nncut0p0_SignalM{m}.root".format(m = mass, var = var))
    t = ROOT.TFile.Open(workdir+"linear_{var}.root".format(var = var), "update")

    for p in process:
      background = f.GetKey("bg_all_{proc}_M{m}".format(m = mass, proc = p))
      data = f.GetKey("data_obs_{proc}_M{m}".format(m = mass, proc = p))
      bg = f.Get(background.GetName())
      da = f.Get(data.GetName())

      bgxBins = bg.GetNbinsX()
      bgyBins = bg.GetNbinsY()
      daxBins = da.GetNbinsX()
      dayBins = da.GetNbinsY()

      if bgxBins != daxBins:
        print "Bad x bins!!!"
      if bgyBins != dayBins:
        print "Bad y bins!!!"

      bgh = ROOT.TH1F("bg_all_{proc}_{var}_M{m}".format(m = mass, proc = p, var = var), "bg_all_{proc}_{var}_M{m}".format(m = mass, proc = p, var = var), bgxBins*bgyBins, 0, 3*bgyBins)
      dah = ROOT.TH1F("data_obs_{proc}_{var}_M{m}".format(m = mass, proc = p, var = var), "data_obs_{proc}_{var}_M{m}".format(m = mass, proc = p, var = var), daxBins*dayBins, 0, 3*dayBins)

      for i in range(bgxBins):
        for j in range(bgyBins):
          bgh.SetBinContent(i*(j+1), bg.Integral(i, i+1, j, j+1))
          dah.SetBinContent(i*(j+1), da.Integral(i, i+1, j, j+1))

      bgh.Write()
      dah.Write()

    total_bgh = ROOT.TH1F("bg_all_ElEl_MuEl_MuMu_{var}_M{m}".format(m = mass, var = var), "bg_all_ElEl_MuEl_MuMu_{var}_M{m}".format(m = mass, var = var), bgxBins*bgyBins, 0, 3*bgyBins)
    total_dah = ROOT.TH1F("data_obs_ElEl_MuEl_MuMu_{var}_M{m}".format(m = mass, var = var), "data_obs_ElEl_MuEl_MuMu_{var}_M{m}".format(m = mass, var = var), daxBins*dayBins, 0, 3*dayBins)
    bgElElh = t.Get(t.GetKey("bg_all_ElEl_{var}_M{m}".format(m = mass, var = var)).GetName())
    bgMuElh = t.Get(t.GetKey("bg_all_MuEl_{var}_M{m}".format(m = mass, var = var)).GetName())
    bgMuMuh = t.Get(t.GetKey("bg_all_MuMu_{var}_M{m}".format(m = mass, var = var)).GetName())
    daElElh = t.Get(t.GetKey("data_obs_ElEl_{var}_M{m}".format(m = mass, var = var)).GetName())
    daMuElh = t.Get(t.GetKey("data_obs_MuEl_{var}_M{m}".format(m = mass, var = var)).GetName())
    daMuMuh = t.Get(t.GetKey("data_obs_MuMu_{var}_M{m}".format(m = mass, var = var)).GetName())

    for i in range(bgxBins):
      for j in range(bgyBins):
        total_bgh.SetBinContent(i*(j+1), bgElElh.GetBinContent(i*(j+1)) + bgMuElh.GetBinContent(i*(j+1)) + bgMuMuh.GetBinContent(i*(j+1)))
        total_dah.SetBinContent(i*(j+1), daElElh.GetBinContent(i*(j+1)) + daMuElh.GetBinContent(i*(j+1)) + daMuMuh.GetBinContent(i*(j+1)))
    total_bgh.Write()
    total_dah.Write()
    t.Close()
