# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: SingleMuPt100_pythia8_cfi --conditions 93X_upgrade2023_realistic_v5 -n 10 --eventcontent FEVTDEBUG --datatier DIGI -s L1 --no_exec --pileup NoPileUp --geometry Extended2023D17 --era Phase2_timing --filein out_digi.root
import FWCore.ParameterSet.Config as cms
import os 

from Configuration.StandardSequences.Eras import eras

process = cms.Process('L1',eras.Phase2_timing)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('out_digi.root'),
    fileNames = cms.untracked.vstring('file:out_GEN_SIM_DIGI.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

### helper for files on dCache/EOS (LPC)
#def useInputDir(process, inputDir, onEOS = True):
#    theInputFiles = []
#    for d in range(len(inputDir)):
#        my_dir = inputDir[d]
#        if not os.path.isdir(my_dir):
#            print "ERROR: This is not a valid directory: ", my_dir
#            if d==len(inputDir)-1:
#                print "ERROR: No input files were selected"
#                exit()
#            continue
#        #print "Proceed to next directory"
#        ls = os.listdir(my_dir)
#        if onEOS:
#            theInputFiles.extend(['file:' + my_dir[:] + x for x in ls if x.endswith('root')])
#        else:
#            ## this works only if you pass the location on pnfs - FIXME for files staring with store/user/...                                                            
#            theInputFiles.extend([my_dir[16:] + x for x in ls if x.endswith('root')])
#
#    process.source.fileNames = cms.untracked.vstring(*theInputFiles)
#inputdir = ['/eos/uscms/store/user/tahuang/SingleMuon_pt2_50_101X_20180610_GEN_SIm_DIGI_93X_upgrade2023_realistic_v5/SingleMuon_pt2_50_101X_20180610_GEN_SIm_DIGI_93X_upgrade2023_realistic_v5/180612_144230/0000/']
#useInputDir(process, inputdir)
#print "inputfiles "
#print process.source.fileNames
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('SingleMuPt100_pythia8_cfi nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DIGI'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('out_L1.root'),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition
process.FEVTDEBUGoutput.outputCommands.append('drop *_simHcal*_*_*')                                                                                       
process.FEVTDEBUGoutput.outputCommands.append('drop *_simEcal*_*_*')
process.FEVTDEBUGoutput.outputCommands.append('drop *_simSiStripDigis_*_*')
process.FEVTDEBUGoutput.outputCommands.append('keep *_emtfStage2Digis_*_*')
process.FEVTDEBUGoutput.outputCommands.append('keep *_simEmtfDigis_*_*')


# For LogTrace to take an effect, compile using
# > scram b -j8 USER_CXXFLAGS="-DEDM_ML_DEBUG"
process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring("log", "debug", "errors"),
    statistics = cms.untracked.vstring("stat"),
    #	untracked vstring categories     = { "lctDigis" }
    #	untracked vstring debugModules   = { "*" }
    #	untracked PSet debugmessages.txt = {
    #	    untracked string threshold = "DEBUG"
    #	    untracked PSet INFO     = {untracked int32 limit = 0}
    #	    untracked PSet DEBUG    = {untracked int32 limit = 0}
    #	    untracked PSet lctDigis = {untracked int32 limit = 10000000}
    #	}
    log = cms.untracked.PSet(
        extension = cms.untracked.string(".txt"),
        lineLength = cms.untracked.int32(132),
        noLineBreaks = cms.untracked.bool(True)
    ),
    debug = cms.untracked.PSet(
        threshold = cms.untracked.string("DEBUG"),
        extension = cms.untracked.string(".txt"),
        lineLength = cms.untracked.int32(132),
        noLineBreaks = cms.untracked.bool(True)
    ),
    errors = cms.untracked.PSet(
        extension = cms.untracked.string(".txt"),
        threshold = cms.untracked.string("ERROR")
    ),
    stat = cms.untracked.PSet(
        extension = cms.untracked.string(".txt"),
        threshold = cms.untracked.string("INFO")
    ),
    debugModules = cms.untracked.vstring("simCscTriggerPrimitiveDigis")
    #debugModules = cms.untracked.vstring("*")
)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '93X_upgrade2023_realistic_v5', '')

# Path and EndPath definitions
process.pL1TkPrimaryVertex = cms.Path(process.L1TkPrimaryVertex)
process.pL1TkElectrons = cms.Path(process.L1TkElectrons)
process.pL1TkJets = cms.Path(process.L1TkJets)
process.pL1TkPhotons = cms.Path(process.L1TkPhotons)
process.pL1TkMuon = cms.Path(process.L1TkMuons)
process.pL1TrkMET = cms.Path(process.L1TkEtMiss)
process.pL1TkTauFromCalo = cms.Path(process.L1TkTauFromCalo)
process.pL1TkHTMissVtx = cms.Path(process.L1TkHTMissVtx)
process.pL1TkIsoElectrons = cms.Path(process.L1TkIsoElectrons)
#process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.simCscTriggerPrimitiveDigis.me11tmbSLHCGEM.verbosity = cms.int32(3)
process.simCscTriggerPrimitiveDigis.me11tmbSLHCGEM.debugMatching = True
process.simCscTriggerPrimitiveDigis.me11tmbSLHCGEM.tmbCrossBxAlgorithm = cms.uint32(1)
process.L1simulation_step = cms.Path(process.simCscTriggerPrimitiveDigis)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.L1simulation_step,process.endjob_step,process.FEVTDEBUGoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
