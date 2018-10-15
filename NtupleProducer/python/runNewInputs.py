import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import FWCore.Utilities.FileUtils as FileUtils

from Configuration.StandardSequences.Eras import eras
process = cms.Process("IN", eras.Phase2_trigger)
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2023_realistic_v1', '')

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff')
process.VertexProducer.l1TracksInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks")

options = VarParsing.VarParsing ('analysis')
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
options.register('winSize',0.33,VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float,"Size of matching window used for vertexing in cm.")
options.register('skip',0,VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int,"No events to skip.")
options.register('job',0,VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int,"No job.")
options.register('outFile', 'ttbar-pu200-937relval-inputs',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,"Output filename prefix.")
options.register('outdir', '.',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,"Output directory.")
options.register('inputFile','',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,"Input filename.")

options.parseArguments()

# matching_window = cms.double(15.05)
matching_window = options.winSize

process.l1pfProducer.vtxRes = matching_window
process.l1pfProducer.AdaptiveCut = cms.bool(True) # whether to treat endcap tracks differently

process.VertexProducer.VertexReconstruction.TP_VertexWidth = matching_window
#process.L1TVertexAnalyzer.VertexReconstruction.TP_VertexWidth = matching_window
process.VertexProducer.VertexReconstruction.VertexResolution = matching_window
#process.L1TVertexAnalyzer.VertexReconstruction.VertexResolution = matching_window

# list = FileUtils.loadListFromFile("937relval-nugun-pu200.txt")
# list = FileUtils.loadListFromFile("937relval-ttbar-pu200.txt")

if (options.inputFile == ''):
    list = FileUtils.loadListFromFile(options.inputFiles)
    readFiles = cms.untracked.vstring(*list)
else:
    readFiles = cms.untracked.vstring(options.inputFile)


process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:/eos/cms/store/relval/CMSSW_9_3_7/RelValSingleTauFlatPt2To100_pythia8/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/58739462-8B2C-E811-A013-0CC47A7C34D0.root'),
    #fileNames = cms.untracked.vstring('file:/eos/cms/store/relval/CMSSW_9_3_7/RelValZMM_14/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v1/10000/96D02123-012D-E811-868C-0CC47A4D7640.root'),
    fileNames = readFiles,
    # fileNames = cms.untracked.vstring('/store/relval/CMSSW_9_3_7/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/7EC7DD7F-782C-E811-B469-0CC47A4D76A0.root'),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    skipEvents = cms.untracked.uint32(options.skip),
    inputCommands = cms.untracked.vstring("keep *", 
        "drop l1tHGCalTowerMapBXVector_hgcalTriggerPrimitiveDigiProducer_towerMap_HLT",
        "drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT",
        "drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT",
        "drop l1tEMTFHit2016s_simEmtfDigis__HLT",
        "drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT",
        "drop l1tEMTFTrack2016s_simEmtfDigis__HLT")

)
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents))
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
  

process.p = cms.Path(
    process.L1TrackletTracks + process.SimL1Emulator
)

process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string("%s/%s-%03i.root" % (options.outdir, options.outFile, options.job)),
        # fileName = cms.untracked.string("inputs.root"),
        outputCommands = cms.untracked.vstring("drop *",
            # --- GEN
            "keep *_genParticles_*_*",
            "keep *_ak4GenJetsNoNu_*_*",
            "keep *_genMetTrue_*_*",
            # --- PF IN
            "keep *_TTTracksFromTracklet_Level1TTTracks_*",
            "keep *_l1EGammaCrystalsProducer_L1EGXtalClusterNoCuts_*",
            "keep *_hgcalTriggerPrimitiveDigiProducer_cluster3D_IN",
            "keep *_hgcalTriggerPrimitiveDigiProducer_cluster2D_IN",
            "keep *_hgcalTriggerPrimitiveDigiProducer_calibratedTriggerCells_IN",
            "keep *_simGmtStage2Digis__*",
            "keep *_simHcalTriggerPrimitiveDigis__*",
            # --- to be tested
            "keep *_hgcalTriggerPrimitiveDigiProducer_tower_IN",
            "keep *_hgcalTriggerPrimitiveDigiProducer_towerMap_IN",
            # --- Stage2 ---
            "keep *_simCaloStage2Digis_*_*",
            # --- TK IN
            "keep *_L1TkCaloHTMissVtx_*_*",
            "keep *_L1TkCaloJets_*_*",
            "keep *_L1TkElectrons_*_*",
            "keep *_L1TkGlbMuons_*_*",
            "keep *_L1TkIsoElectrons_*_*",
            "keep *_L1TkMuons_*_*",
            "keep *_L1TkPhotons_*_*",
            "keep *_L1TkPrimaryVertex_*_*",
            "keep *_L1TkTauFromCalo_*_*",
            "keep *_L1TrackerEtMiss_*_*",
            "keep *_L1TrackerHTMiss_*_*",
            "keep *_L1TrackerJets_*_*",
            "keep *_VertexProducer_*_*",
            # --- PF OUT
            "keep l1tPFClusters_*_*_*",
            "keep l1tPFTracks_*_*_*",
            "keep l1tPFCandidates_*_*_*",
        ),
        compressionAlgorithm = cms.untracked.string('LZMA'),
        compressionLevel = cms.untracked.int32(4),
        dropMetaData = cms.untracked.string('ALL'),
        fastCloning = cms.untracked.bool(False),
        overrideInputFileSplitLevels = cms.untracked.bool(True),
        eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
)
process.e = cms.EndPath(process.out)

process.schedule = cms.Schedule([process.p,process.e])
