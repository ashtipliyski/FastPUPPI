import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from math import sqrt

process = cms.Process("RESP", eras.Phase2C4_trigger)

process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/eos/cms/store/cmst3/user/gpetrucc/l1tr/105X/NewInputs104X/010319/TTbar_PU200/inputs104X_TTbar_PU200_job1.root'),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    skipBadFiles = cms.untracked.bool(True)
)

process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D35_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMET.Configuration.GenMETParticles_cff')
process.load('RecoMET.METProducers.genMetTrue_cfi')

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
from RecoMET.METProducers.PFMET_cfi import pfMet

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '') 

process.load("L1Trigger.Phase2L1ParticleFlow.l1ParticleFlow_cff")
process.load("L1Trigger.Phase2L1ParticleFlow.l1ParticleFlow_old_cff")

process.extraPFStuff = cms.Task()

process.runPF = cms.Sequence( 
    process.l1ParticleFlow 
)

process.centralGen = cms.EDFilter("CandPtrSelector", src = cms.InputTag("genParticlesForMETAllVisible"), cut = cms.string("abs(eta) < 2.4"))
process.barrelGen = cms.EDFilter("CandPtrSelector", src = cms.InputTag("genParticlesForMETAllVisible"), cut = cms.string("abs(eta) < 1.5"))
process.genMetCentralTrue = process.genMetTrue.clone(src = cms.InputTag("centralGen"))
process.genMetBarrelTrue = process.genMetTrue.clone(src = cms.InputTag("barrelGen"))
process.extraPFStuff.add(
    process.genParticlesForMETAllVisible,
    process.centralGen,
    process.barrelGen,
    process.genMetCentralTrue,
    process.genMetBarrelTrue
)

def monitorPerf(label, tag, makeResp=True, makeRespSplit=True, makeJets=True, makeMET=True, makeCentralMET=True, makeBarrelMET=True):
    def _add(name, what):
        setattr(process, name, what)
        process.extraPFStuff.add(what)
    if type(tag) != str and len(tag) > 1:
        _add('merged'+label, cms.EDProducer("L1TPFCandMerger", src = cms.VInputTag(cms.InputTag(x) for x in tag)))
        tag = 'merged'+label
    if makeResp:
        setattr(process.ntuple.objects, label, cms.VInputTag(cms.InputTag(tag)))
        if makeRespSplit:
            setattr(process.ntuple.objects, label+"Charged", cms.VInputTag(cms.InputTag(tag)))
            setattr(process.ntuple.objects, label+"Charged_sel", cms.string("charge != 0"))
            setattr(process.ntuple.objects, label+"Photon",  cms.VInputTag(cms.InputTag(tag)))
            setattr(process.ntuple.objects, label+"Photon_sel", cms.string("pdgId == 22"))
            setattr(process.ntuple.objects, label+"Neutral",  cms.VInputTag(cms.InputTag(tag)))
            setattr(process.ntuple.objects, label+"Neutral_sel", cms.string("charge == 0"))
            setattr(process.ntuple.objects, label+"NeutralHad",  cms.VInputTag(cms.InputTag(tag)))
            setattr(process.ntuple.objects, label+"NeutralHad_sel", cms.string("charge == 0 && pdgId != 22"))
    if makeJets:
        _add('ak4'+label, ak4PFJets.clone(src = tag, doAreaFastjet = False))
        setattr(process.l1pfjetTable.jets, label, cms.InputTag('ak4'+label))
    if makeMET:
        _add('met'+label, pfMet.clone(src = tag, calculateSignificance = False))
        setattr(process.l1pfmetTable.mets, label, cms.InputTag('met'+label))
        if makeCentralMET:
            _add('central'+label, cms.EDFilter("CandPtrSelector", src = cms.InputTag(tag), cut = cms.string("abs(eta) < 2.4")))
            _add('met'+label+'Central', pfMet.clone(src = 'central'+label, calculateSignificance = False))
            setattr(process.l1pfmetCentralTable.mets, label, cms.InputTag('met'+label+'Central'))
        if makeBarrelMET:
            _add('barrel'+label, cms.EDFilter("CandPtrSelector", src = cms.InputTag(tag), cut = cms.string("abs(eta) < 1.5")))
            _add('met'+label+'Barrel', pfMet.clone(src = 'barrel'+label, calculateSignificance = False))
            setattr(process.l1pfmetBarrelTable.mets, label, cms.InputTag('met'+label+'Barrel'))


class GatherAllModulesVisitor(object):
    """Visitor that travels within a cms.Sequence, and returns a list of objects of type gatheredInance(e.g. modules) that have it"""
    def __init__(self, gatheredInstance=cms._Module):
        self._modules = []
        self._gatheredInstance= gatheredInstance
    def enter(self,visitee):
        if isinstance(visitee,self._gatheredInstance):
            self._modules.append(visitee)
    def leave(self,visitee):
        pass
    def modules(self):
        return self._modules

def listModules(sequence):
    visitor = GatherAllModulesVisitor(gatheredInstance=cms._Module)
    sequence.visit(visitor)
    return visitor.modules()

def listSequences(sequence):
    visitor = GatherAllModulesVisitor(gatheredInstance=cms.Sequence)
    sequence.visit(visitor)
    return visitor.modules()

def __labelsInSequence(process, sequenceLabel, postfix="", keepPostFix=False):
    position = -len(postfix)
    if keepPostFix: 
        position = None

    result = [ m.label()[:position] for m in listModules( getattr(process,sequenceLabel+postfix))]
    result.extend([ m.label()[:position] for m in listSequences( getattr(process,sequenceLabel+postfix))]  )
    if postfix == "":
        result = [ m.label() for m in listModules( getattr(process,sequenceLabel+postfix))]
        result.extend([ m.label() for m in listSequences( getattr(process,sequenceLabel+postfix))]  )
    return result



process.ntuple = cms.EDAnalyzer("ResponseNTuplizer",
    genJets = cms.InputTag("ak4GenJetsNoNu"),
    genParticles = cms.InputTag("genParticles"),
    isParticleGun = cms.bool(False),
    doRandom = cms.bool(False),
    objects = cms.PSet(
        # -- inputs and PF --
        RawTK  = cms.VInputTag('pfTracksFromL1Tracks',),
        # outputs
    ),
    copyUInts = cms.VInputTag(),
    copyFloats = cms.VInputTag(),
)
process.extraPFStuff.add(process.pfTracksFromL1Tracks)

process.l1pfjetTable = cms.EDProducer("L1PFJetTableProducer",
    gen = cms.InputTag("ak4GenJetsNoNu"),
    commonSel = cms.string("pt > 5 && abs(eta) < 5.0"),
    drMax = cms.double(0.2),
    minRecoPtOverGenPt = cms.double(0.1),
    jets = cms.PSet(
        Gen = cms.InputTag("ak4GenJetsNoNu"),
        Gen_sel = cms.string("pt > 15"),
    ),
    moreVariables = cms.PSet(
    ),
)

process.l1pfmetTable = cms.EDProducer("L1PFMetTableProducer",
    genMet = cms.InputTag("genMetTrue"), 
    flavour = cms.string(""),
    mets = cms.PSet(
    ),
)

process.l1ParticleFlow1 = cms.Sequence()
for l_ in  __labelsInSequence(process,"l1ParticleFlow"):
    if (not hasattr(getattr(process,l_),"_seq")):
        newl = getattr(process,l_).clone()
        setattr(process,l_+"1",newl)
        process.l1ParticleFlow1 += getattr(process,l_+"1")

process.l1ParticleFlow1.replace("l1pfProducerBarrel","l1pfProducerBarrel1")
process.l1ParticleFlow1.replace("l1pfProducerHGCal","l1pfProducerHGCal1")
process.l1ParticleFlow1.replace("l1pfProducerHF","l1pfProducerHF1")
process.l1pfProducerBarrel1.vtxNum = cms.untracked.uint32(1)
process.l1pfProducerHGCal1.vtxNum = cms.untracked.uint32(1)
process.l1pfProducerHF1.vtxNum = cms.untracked.uint32(1)
process.l1pfCandidates1.pfProducers = cms.VInputTag(
        cms.InputTag("l1pfProducerBarrel1"),
        cms.InputTag("l1pfProducerHGCal1"),
        cms.InputTag("l1pfProducerHF1")
    )
process.l1ParticleFlow1 = cms.Sequence()
for l_ in  __labelsInSequence(process,"l1ParticleFlow"):
    if (not hasattr(getattr(process,l_),"_seq")):
        newl = getattr(process,l_).clone()
        setattr(process,l_+"1",newl)
        process.l1ParticleFlow1 += getattr(process,l_+"1")

process.l1ParticleFlow1.replace("l1pfProducerBarrel","l1pfProducerBarrel1")
process.l1ParticleFlow1.replace("l1pfProducerHGCal","l1pfProducerHGCal1")
process.l1ParticleFlow1.replace("l1pfProducerHF","l1pfProducerHF1")
process.l1pfProducerBarrel1.vtxNum = cms.untracked.uint32(1)
process.l1pfProducerHGCal1.vtxNum = cms.untracked.uint32(1)
process.l1pfProducerHF1.vtxNum = cms.untracked.uint32(1)
process.l1pfCandidates1.pfProducers = cms.VInputTag(
        cms.InputTag("l1pfProducerBarrel1"),
        cms.InputTag("l1pfProducerHGCal1"),
        cms.InputTag("l1pfProducerHF1")
    )


process.l1ParticleFlow2 = cms.Sequence()
for l_ in  __labelsInSequence(process,"l1ParticleFlow"):
    if (not hasattr(getattr(process,l_),"_seq")):
        newl = getattr(process,l_).clone()
        setattr(process,l_+"2",newl)
        process.l1ParticleFlow2 += getattr(process,l_+"2")

process.l1ParticleFlow2.replace("l1pfProducerBarrel","l1pfProducerBarrel2")
process.l1ParticleFlow2.replace("l1pfProducerHGCal","l1pfProducerHGCal2")
process.l1ParticleFlow2.replace("l1pfProducerHF","l1pfProducerHF2")
process.l1pfProducerBarrel2.vtxNum = cms.untracked.uint32(2)
process.l1pfProducerHGCal2.vtxNum = cms.untracked.uint32(2)
process.l1pfProducerHF2.vtxNum = cms.untracked.uint32(2)
process.l1pfCandidates2.pfProducers = cms.VInputTag(
        cms.InputTag("l1pfProducerBarrel2"),
        cms.InputTag("l1pfProducerHGCal2"),
        cms.InputTag("l1pfProducerHF2")
    )


process.l1ParticleFlow3 = cms.Sequence()
for l_ in  __labelsInSequence(process,"l1ParticleFlow"):
    if (not hasattr(getattr(process,l_),"_seq")):
        newl = getattr(process,l_).clone()
        setattr(process,l_+"3",newl)
        process.l1ParticleFlow3 += getattr(process,l_+"3")

process.l1ParticleFlow3.replace("l1pfProducerBarrel","l1pfProducerBarrel3")
process.l1ParticleFlow3.replace("l1pfProducerHGCal","l1pfProducerHGCal3")
process.l1ParticleFlow3.replace("l1pfProducerHF","l1pfProducerHF3")
process.l1pfProducerBarrel3.vtxNum = cms.untracked.uint32(3)
process.l1pfProducerHGCal3.vtxNum = cms.untracked.uint32(3)
process.l1pfProducerHF3.vtxNum = cms.untracked.uint32(3)
process.l1pfCandidates3.pfProducers = cms.VInputTag(
        cms.InputTag("l1pfProducerBarrel3"),
        cms.InputTag("l1pfProducerHGCal3"),
        cms.InputTag("l1pfProducerHF3")
    )

process.l1ParticleFlow4 = cms.Sequence()
for l_ in  __labelsInSequence(process,"l1ParticleFlow"):
    if (not hasattr(getattr(process,l_),"_seq")):
        newl = getattr(process,l_).clone()
        setattr(process,l_+"4",newl)
        process.l1ParticleFlow4 += getattr(process,l_+"4")

process.l1ParticleFlow4.replace("l1pfProducerBarrel","l1pfProducerBarrel4")
process.l1ParticleFlow4.replace("l1pfProducerHGCal","l1pfProducerHGCal4")
process.l1ParticleFlow4.replace("l1pfProducerHF","l1pfProducerHF4")
process.l1pfProducerBarrel4.vtxNum = cms.untracked.uint32(4)
process.l1pfProducerHGCal4.vtxNum = cms.untracked.uint32(4)
process.l1pfProducerHF4.vtxNum = cms.untracked.uint32(4)
process.l1pfCandidates4.pfProducers = cms.VInputTag(
        cms.InputTag("l1pfProducerBarrel4"),
        cms.InputTag("l1pfProducerHGCal4"),
        cms.InputTag("l1pfProducerHF4")
    )



process.multweights = cms.EDProducer("L1MultPuppiWeightProducer",
    pf = cms.InputTag("l1pfCandidates:PF"),
    pup = cms.InputTag("l1pfCandidates:Puppi"),
    pup1 = cms.InputTag("l1pfCandidates1:Puppi"),
    pup2 = cms.InputTag("l1pfCandidates2:Puppi"),
    pup3 = cms.InputTag("l1pfCandidates3:Puppi"),
    pup4 = cms.InputTag("l1pfCandidates4:Puppi"),
    pup5 = cms.InputTag("l1pfCandidates5:Puppi"),
    top_m = cms.uint32(1)
)


process.l1pfmetCentralTable = process.l1pfmetTable.clone(genMet = "genMetCentralTrue", flavour = "Central")
process.l1pfmetBarrelTable  = process.l1pfmetTable.clone(genMet = "genMetBarrelTrue", flavour = "Barrel")

monitorPerf("L1Calo", "l1pfCandidates:Calo", makeRespSplit = False)
monitorPerf("L1TK", "l1pfCandidates:TK", makeRespSplit = False, makeJets=False, makeMET=False)
monitorPerf("L1TKV", "l1pfCandidates:TKVtx", makeRespSplit = False, makeJets=False, makeMET=False)
monitorPerf("L1PF", "l1pfCandidates:PF")
# monitorPerf("L1Puppi", "multweights:MultiPuppi")
monitorPerf("L1Puppi", "l1pfCandidates:Puppi")

process.runPF.associate(process.extraPFStuff)
process.p = cms.Path(
        process.runPF + 
        process.ntuple + #
        process.l1ParticleFlow1 +
        process.l1ParticleFlow2 +
        process.l1ParticleFlow3 +
        process.l1ParticleFlow4 +
        process.multweights +
        process.l1pfjetTable + 
        process.l1pfmetTable + process.l1pfmetCentralTable + process.l1pfmetBarrelTable
        )
process.TFileService = cms.Service("TFileService", fileName = cms.string("perfTuple.root"))

process.outnano = cms.OutputModule("NanoAODOutputModule",
    fileName = cms.untracked.string("perfNano.root"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    outputCommands = cms.untracked.vstring("drop *", "keep nanoaodFlatTable_*Table_*_*"),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("ZLIB"),
)
process.end = cms.EndPath(process.outnano)

# Below for more debugging
if True:
    process.genInAcceptance = cms.EDFilter("GenParticleSelector",
        src = cms.InputTag("genParticles"),
        cut = cms.string("status == 1 && (abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16) && "+
                         "(abs(eta) < 2.5 && pt > 2 && charge != 0 || "+
                         "abs(pdgId) == 22 && pt > 1 || "+
                         "charge == 0 && pt > 1 || "+
                         "charge != 0 && abs(eta) > 2.5 && pt > 2) ") # tracks below pT 2 bend by more than 0.4,
    )
    process.ntuple.objects.GenAcc = cms.VInputTag(cms.InputTag("genInAcceptance"))
    process.ntuple.objects.ChGenAcc = cms.VInputTag(cms.InputTag("genInAcceptance"))
    process.ntuple.objects.ChGenAcc_sel = cms.string("(abs(eta) < 2.5 && pt > 2 && charge != 0)")
    process.ntuple.objects.PhGenAcc = cms.VInputTag(cms.InputTag("genInAcceptance"))
    process.ntuple.objects.PhGenAcc_sel = cms.string("pdgId == 22")
    process.extraPFStuff.add(process.genInAcceptance)
if False: # test also PF leptons
    process.ntuple.objects.L1PFMuon = cms.VInputTag("l1pfCandidates:PF",)
    process.ntuple.objects.L1PFMuon_sel = cms.string("abs(pdgId) == 13")
    process.ntuple.objects.L1PFElectron = cms.VInputTag("l1pfCandidates:PF",)
    process.ntuple.objects.L1PFElectron_sel = cms.string("abs(pdgId) == 11")
def respOnly():
    process.p.remove(process.l1pfjetTable)
    process.p.remove(process.l1pfmetTable)
    process.p.remove(process.l1pfmetCentralTable)
    process.p.remove(process.l1pfmetBarrelTable)
    process.end.remove(process.outnano)
def addCHS():
    process.l1PuppiCharged = cms.EDFilter("L1TPFCandSelector",
        src = cms.InputTag("l1pfCandidates:Puppi"),
        cut = cms.string("charge != 0"))
    process.l1PFNeutral = cms.EDFilter("L1TPFCandSelector",
        src = cms.InputTag("l1pfCandidates:PF"),
        cut = cms.string("charge == 0"))
    process.extraPFStuff.add(process.l1PuppiCharged, process.l1PFNeutral)
    monitorPerf("L1CHS", [ "l1PuppiCharged", "l1PFNeutral" ], makeRespSplit = False)
def addOld():
    process.extraPFStuff.add(
        process.pfClustersFromHGC3DClustersEM, 
        process.pfClustersFromCombinedCalo, 
        process.l1pfProducerOld,
        process.l1OldPuppiForMET
    )
    monitorPerf("L1OldCalo", "l1pfProducerOld:Calo", makeRespSplit = False)
    monitorPerf("L1OldPF", "l1pfProducerOld:PF")
    monitorPerf("L1OldPuppi", "l1pfProducerOld:Puppi")
    monitorPerf("L1OldPuppiForMET", "l1OldPuppiForMET")
def addPuppiOld():
    process.l1pfProducerBarrelPuppiOld = process.l1pfProducerBarrel.clone(puAlgo = "Puppi")
    process.l1pfProducerBarrelPuppiOld.puppiEtaCuts       = cms.vdouble(1.5) 
    process.l1pfProducerBarrelPuppiOld.puppiPtCuts        = cms.vdouble(0.0)
    process.l1pfProducerBarrelPuppiOld.puppiPtCutsPhotons = cms.vdouble(0.0)
    #
    process.l1pfProducerHGCalPuppiOld = process.l1pfProducerHGCal.clone(puAlgo = "Puppi")
    process.l1pfProducerHGCalPuppiOld.puppiEtaCuts       = cms.vdouble(2.5, 2.85,  3.0) 
    process.l1pfProducerHGCalPuppiOld.puppiPtCuts        = cms.vdouble( 20,  40,  9999)
    process.l1pfProducerHGCalPuppiOld.puppiPtCutsPhotons = cms.vdouble( 20,  40,  9999)
    #
    process.l1pfProducerHFPuppiOld = process.l1pfProducerHF.clone(puAlgo = "Puppi")
    process.l1pfProducerHFPuppiOld.puppiEtaCuts       = cms.vdouble(5.5) 
    process.l1pfProducerHFPuppiOld.puppiPtCuts        = cms.vdouble( 30)
    process.l1pfProducerHFPuppiOld.puppiPtCutsPhotons = cms.vdouble( 30)
    #
    process.l1pfCandidatesPuppiOld = process.l1pfCandidates.clone(
            pfProducers = ["l1pfProducerBarrelPuppiOld", "l1pfProducerHGCalPuppiOld", "l1pfProducerHFPuppiOld"], 
            labelsToMerge = [ "Puppi" ])
    process.extraPFStuff.add(process.l1pfProducerBarrelPuppiOld, process.l1pfProducerHGCalPuppiOld, process.l1pfProducerHFPuppiOld, process.l1pfCandidatesPuppiOld)
    monitorPerf("L1PuppiOld", "l1pfCandidatesPuppiOld:Puppi")
def addTKs():
    process.l1tkv5Stubs = cms.EDFilter("L1TPFCandSelector", src = cms.InputTag("l1pfCandidates:TKVtx"), cut = cms.string("pfTrack.nStubs >= 5"))
    process.l1tkv6Stubs = cms.EDFilter("L1TPFCandSelector", src = cms.InputTag("l1pfCandidates:TKVtx"), cut = cms.string("pfTrack.nStubs >= 6"))
    process.extraPFStuff.add(process.l1tkv5Stubs, process.l1tkv6Stubs)
    monitorPerf("L1TKV5", "l1tkv5Stubs", makeRespSplit = False)
    monitorPerf("L1TKV6", "l1tkv6Stubs", makeRespSplit = False)
    monitorPerf("L1TK", "l1pfCandidates:TK", makeRespSplit = False)
    monitorPerf("L1TKV", "l1pfCandidates:TKVtx", makeRespSplit = False)
def addCalib():
    process.pfClustersFromL1EGClustersRaw    = process.pfClustersFromL1EGClusters.clone(corrector = "")
    process.pfClustersFromHGC3DClustersRaw   = process.pfClustersFromHGC3DClusters.clone(corrector = "")
    process.pfClustersFromHGC3DClustersEMRaw = process.pfClustersFromHGC3DClustersRaw.clone(emOnly = True, etMin = 0.)
    process.extraPFStuff.add(
            process.pfClustersFromL1EGClustersRaw, 
            process.pfClustersFromHGC3DClustersRaw, 
            process.pfClustersFromHGC3DClustersEMRaw)
    process.ntuple.objects.L1RawBarrelEcal   = cms.VInputTag('pfClustersFromL1EGClustersRaw' )
    process.ntuple.objects.L1RawBarrelCalo   = cms.VInputTag('pfClustersFromCombinedCaloHCal:uncalibrated')
    process.ntuple.objects.L1RawBarrelCaloEM = cms.VInputTag('pfClustersFromCombinedCaloHCal:emUncalibrated')
    process.ntuple.objects.L1RawHGCal   = cms.VInputTag('pfClustersFromHGC3DClustersRaw')
    process.ntuple.objects.L1RawHGCalEM = cms.VInputTag('pfClustersFromHGC3DClustersEMRaw')
    process.ntuple.objects.L1RawHFCalo  = cms.VInputTag('pfClustersFromCombinedCaloHF:uncalibrated')
    process.ntuple.objects.L1BarrelEcal = cms.VInputTag('pfClustersFromL1EGClusters' )
    process.ntuple.objects.L1BarrelCalo = cms.VInputTag('pfClustersFromCombinedCaloHCal:calibrated')
    process.ntuple.objects.L1HGCal   = cms.VInputTag('pfClustersFromHGC3DClusters')
    process.ntuple.objects.L1HFCalo  = cms.VInputTag('pfClustersFromCombinedCaloHF:calibrated')
    if hasattr(process.ntuple.objects, 'L1OldPF'):
        process.ntuple.objects.L1HGCalEM = cms.VInputTag('pfClustersFromHGC3DClustersEM', )
        process.ntuple.objects.L1OldRawCalo = cms.VInputTag('pfClustersFromCombinedCalo:uncalibrated')
        process.ntuple.objects.L1OldRawCaloEM = cms.VInputTag('pfClustersFromCombinedCalo:emUncalibrated')
        process.ntuple.objects.L1OldRawEcal = cms.VInputTag('pfClustersFromL1EGClustersRaw', 'pfClustersFromHGC3DClustersEMRaw')
        process.ntuple.objects.L1OldEcal = cms.VInputTag(cms.InputTag('l1pfProducerOld','EmCalo'))

def goGun(calib=1):
    process.ntuple.isParticleGun = True
    respOnly()
    if calib: 
        addCalib()
def goMT(nthreads=2):
    process.options.numberOfThreads = cms.untracked.uint32(nthreads)
    process.options.numberOfStreams = cms.untracked.uint32(0)
def goOld():
    process.pfClustersFromL1EGClusters.corrector  = "L1Trigger/Phase2L1ParticleFlow/data/emcorr_barrel_93X.root"
    process.pfClustersFromHGC3DClustersEM.corrector =  "L1Trigger/Phase2L1ParticleFlow/data/emcorr_hgc_old3d_93X.root"
    process.pfClustersFromCombinedCalo.hadCorrector =  "L1Trigger/Phase2L1ParticleFlow/data/hadcorr_93X.root"
    process.pfClustersFromHGC3DClusters.corrector =  "L1Trigger/Phase2L1ParticleFlow/data/hadcorr_HGCal3D_STC_93X.root"
    process.pfClustersFromCombinedCaloHCal.hadCorrector =  "L1Trigger/Phase2L1ParticleFlow/data/hadcorr_barrel_93X.root"
    process.pfClustersFromCombinedCaloHF.hadCorrector =  "L1Trigger/Phase2L1ParticleFlow/data/hfcorr_93X.root"

