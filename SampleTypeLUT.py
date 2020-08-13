
sampletypedict = {
		"TT":13,
		"TTbar":13,
		"DYM10to50":14,
		"DYToLL0J":15,
		"DYToLL1J":16,
		"DYToLL2J":17,
		"sT_top":18,
		"sT_antitop":19,
		}
resonant_allmasses = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900]
for i,m in enumerate(resonant_allmasses):
	sampletypedict["radion_M%d"%m] = 100+i
	sampletypedict["graviton_M%d"%m] = 100+i+len(resonant_allmasses)



            #'workingdir' : '/Users/taohuang/Documents/DiHiggs/20170530/20170905_Louvain_10k/HH_Keras/HHTools/mvaTraining/hh_resonant_trained_models_kinematicwithMT_20171102/2017-11-08_90epochs/',
ModelLUT = { 
        #'MTonly':{
	#    'workingdir' : '/Users/taohuang/Documents/DiHiggs/20170530/20170905_Louvain_10k/HH_Keras/HHTools/mvaTraining/hh_resonant_trained_models_kinematicwithMT_20171110_Louvain/2017-11-10_100epochs/',
        #    'legend': 'Kinematic + MT',
	#    'inputs': ['jj_pt', 'll_pt', 'll_M', 'll_DR_l_l', 'jj_DR_j_j', 'llmetjj_DPhi_ll_jj', 'llmetjj_minDR_l_j', 'llmetjj_MTformula', 'isSF'],
	#	},
	'MTonly':{
	    'workingdir' : '/Users/taohuang/Documents/DiHiggs/HH_Keras/HHTools/mvaTraining/hh_resonant_trained_models_kinematicwithMTonly_20180206/2018-02-06_100epochs/',
            'legend': 'Kinematic + MT',#2018
	    'inputs': ['jj_pt', 'll_pt', 'll_M', 'll_DR_l_l', 'jj_DR_j_j', 'llmetjj_DPhi_ll_jj', 'llmetjj_minDR_l_j', 'llmetjj_MTformula', 'isSF'],
		},
	'MT2only':{
	    'workingdir' : '/Users/taohuang/Documents/DiHiggs/20170530/20170905_Louvain_10k/HH_Keras/HHTools/mvaTraining/hh_resonant_trained_models_kinematicWithMT2withoutMT/2017-10-25_90epochs/',
            'legend': 'Kinematic + MT2',
	    'inputs': ['jj_pt', 'll_pt', 'll_M', 'll_DR_l_l', 'jj_DR_j_j', 'llmetjj_DPhi_ll_jj', 'llmetjj_minDR_l_j', 'mt2', 'isSF'],
		},
        #'MTandMT2':{
	#    'workingdir' : '/Users/taohuang/Documents/DiHiggs/20170530/20170905_Louvain_10k/HH_Keras/HHTools/mvaTraining/hh_resonant_trained_models_kinematicWithMT2andMT/2017-10-25_90epochs/',
        #    'legend': 'Kinematic + MT + MT2',
	#    'inputs': ['jj_pt', 'll_pt', 'll_M', 'll_DR_l_l', 'jj_DR_j_j', 'llmetjj_DPhi_ll_jj', 'llmetjj_minDR_l_j', 'llmetjj_MTformula','mt2', 'isSF'],
        #	},
	'MTandMT2':{
	    'workingdir' : '/Users/taohuang/Documents/DiHiggs/HH_Keras/HHTools/mvaTraining/hh_resonant_trained_models_kinematicwithMTandMT2_20180206/2018-02-07_100epochs/',
            'legend': 'Kinematic + MT + MT2',#2018
	    'inputs': ['jj_pt', 'll_pt', 'll_M', 'll_DR_l_l', 'jj_DR_j_j', 'llmetjj_DPhi_ll_jj', 'llmetjj_minDR_l_j', 'llmetjj_MTformula','mt2','isSF'],
		},
	'MTandMJJ':{
	    'workingdir' : '/Users/taohuang/Documents/DiHiggs/20170530/20170905_Louvain_10k/HH_Keras/HHTools/mvaTraining/hh_resonant_trained_models_kinematicwithMTandMJJ_20171102/2017-11-08_90epochs/',
            'legend': 'Kinematic + MT + MJJ',
	    'inputs': ['jj_pt', 'll_pt', 'll_M', 'll_DR_l_l', 'jj_DR_j_j', 'llmetjj_DPhi_ll_jj', 'llmetjj_minDR_l_j', 'llmetjj_MTformula', 'isSF', 'jj_M'],
	    },
        #'MTandMT2_MJJ':{ 
	#    'workingdir' : '/Users/taohuang/Documents/DiHiggs/20170530/20170905_Louvain_10k/HH_Keras/HHTools/mvaTraining/hh_resonant_trained_models_kinematicWithMT2andMT_Mjj/2017-10-26_90epochs/',
	#    'legend': 'Kinematic+MT+MT2+Mjj',
	#    'inputs': ['jj_pt', 'll_pt', 'll_M', 'll_DR_l_l', 'jj_DR_j_j', 'llmetjj_DPhi_ll_jj', 'llmetjj_minDR_l_j', 'llmetjj_MTformula', 'mt2', 'isSF', 'jj_M'],
        #    },
	'MTandMT2_MJJ':{ 
	    'workingdir' : '/Users/taohuang/Documents/DiHiggs/HH_Keras/HHTools/mvaTraining/hh_resonant_trained_models_kinematicwithMTandMT2Mjj_20180206/2018-02-07_100epochs/',
	    'legend': 'Kinematic+MT+MT2+Mjj',#2018
	    'inputs': ['jj_pt', 'll_pt', 'll_M', 'll_DR_l_l', 'jj_DR_j_j', 'llmetjj_DPhi_ll_jj', 'llmetjj_minDR_l_j', 'llmetjj_MTformula', 'mt2', 'isSF', 'jj_M'],
	    },
	'MTandMT2_HME':{ 
	    'workingdir' : '/Users/taohuang/Documents/DiHiggs/20170530/20170905_Louvain_10k/HH_Keras/HHTools/mvaTraining/hh_resonant_trained_models_kinematicWithMT2andMT_HMERecoMass/2017-10-26_90epochs/',
	    'legend': 'Kinematic+MT+MT2+HME',
	    'inputs': ['jj_pt', 'll_pt', 'll_M', 'll_DR_l_l', 'jj_DR_j_j', 'llmetjj_DPhi_ll_jj', 'llmetjj_minDR_l_j', 'llmetjj_MTformula', 'mt2', 'isSF', 'hme_h2mass_reco'],
	    },
        #'MTandMT2_HMEMJJ':{ 
	#    'workingdir' : '/Users/taohuang/Documents/DiHiggs/20170530/20170905_Louvain_10k/HH_Keras/HHTools/mvaTraining/hh_resonant_trained_models_kinematicWithMT2andMT_HMERecoMassAndMJJ/2017-10-26_90epochs/',
	#    'legend': 'Kinematic+MT+MT2+Mjj+HME',
	#    'inputs': ['jj_pt', 'll_pt', 'll_M', 'll_DR_l_l', 'jj_DR_j_j', 'llmetjj_DPhi_ll_jj', 'llmetjj_minDR_l_j', 'llmetjj_MTformula', 'mt2', 'isSF', 'jj_M', 'hme_h2mass_reco'],
	#    },
	'MTandMT2_HMEMJJ':{ 
	    'workingdir' : '/Users/taohuang/Documents/DiHiggs/HH_Keras/HHTools/mvaTraining/hh_resonant_trained_models_kinematicwithMTandMT2MjjHME_20180206/2018-02-07_100epochs/',
	    'legend': 'Kinematic+MT+MT2+Mjj+HME',#2018
	    'inputs': ['jj_pt', 'll_pt', 'll_M', 'll_DR_l_l', 'jj_DR_j_j', 'llmetjj_DPhi_ll_jj', 'llmetjj_minDR_l_j', 'llmetjj_MTformula', 'mt2', 'isSF', 'jj_M', 'hme_h2mass_reco'],
	    }

	}

def getSignalMass(itype):
	#if type(sampletye) == type.string
	#itype = sampletypedict[sampletype]

	x = itype-100
	return resonant_allmasses[x%len(resonant_allmasses)]
