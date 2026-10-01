# End-to-end Rivet equivalence test that generates Pythia8 ttbar events
# from scratch (no input file), runs BOTH HepMC3 conversion paths in
# parallel, and feeds each to its own RivetAnalyzer so the two sets of
# Rivet plots can be compared.
#
# The two paths, run on the SAME generated events:
#
#   (1) Standard path (what Rivet used in CMSSW so far):
#         Pythia8 -> GenParticleProducer -> GenParticles2HepMCConverter
#                                            -> rivetAnalyzer (-> std.yoda)
#
#   (2) Truth-graph path (the new "Truth Graph" feature):
#         Pythia8 -> EmptySimContainersProducer  (so TruthGraphProducer's
#                    evt.get(SimTracks/SimVertices) finds the branches; with
#                    empty containers the producer builds a GEN-only graph)
#                 -> TruthGraphProducer -> TruthLogicalGraphProducer
#                 -> TruthGraph2HepMCConverter
#                 -> rivetAnalyzerTruthGraph (-> truth.yoda)
#
# The two RivetAnalyzers are configured identically (same analysis names,
# same output options) and run on the same generated events, so the two YODA
# files they write are directly comparable. Use the compareYodaFiles.py helper
# to verify the plots match.
#
# Usage:
#   cmsRun test/runRivetEquivalence_cfg.py [-n NEVTS] [--analyses A,B,...]
#                                        [--std-out FILE] [--truth-out FILE]
#
# Defaults: 100 events, analyses MC_TTBAR and CMS_2018_I1662081, outputs
# mc_rivet_std.yoda / mc_rivet_truth.yoda in the current directory.

import os
import argparse
import FWCore.ParameterSet.Config as cms

_parser = argparse.ArgumentParser(description="Rivet equivalence test (standard vs truth-graph converter)")
_parser.add_argument("-n", "--maxEvents", type=int, default=100,
                    help="number of Pythia8 events to generate (default: %(default)s)")
_parser.add_argument("--analyses", default="MC_TTBAR,MC_FSPARTICLES,MC_PARTONICTOPS,MC_HFDECAYS,MC_JETS",
                    help="comma-separated Rivet analysis names (default: %(default)s)")
_parser.add_argument("--std-out", default="mc_rivet_std.yoda",
                    help="YODA output file of the standard path (default: %(default)s)")
_parser.add_argument("--truth-out", default="mc_rivet_truth.yoda",
                    help="YODA output file of the truth-graph path (default: %(default)s)")
_parser.add_argument("--comEnergy", type=float, default=14000.,
                    help="Pythia8 centre-of-mass energy in GeV (default: %(default)s)")
_args, _unknown = _parser.parse_known_args()

_analysisNames = [name for name in _args.analyses.split(",") if name]

process = cms.Process("RIVETEQ")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.SimulationRandomNumberGeneratorSeeds_cff")
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(_args.maxEvents))
process.options = cms.untracked.PSet(
    allowUnscheduled=cms.untracked.bool(True),
    wantSummary=cms.untracked.bool(True),
)
process.options.numberOfConcurrentLuminosityBlocks = 1

# --- Common Pythia8 ttbar generator (single source for both paths) --------
# The Pythia8GeneratorFilter produces, per event:
#   - edm::HepMCProduct:unsmeared     (HepMC2 record, used by TruthGraphProducer)
#   - edm::HepMC3Product:unsmeared     (HepMC3 record, used by TruthLogicalGraphProducer)
#   - GenEventInfoProduct             (signal process id, weights, PDF info, ...)
# Both conversion paths read from this same generator, so they see the
# exact same generated events.
process.generator = cms.EDFilter("Pythia8GeneratorFilter",
    maxEventsToPrint=cms.untracked.int32(0),
    pythiaPylistVerbosity=cms.untracked.int32(0),
    filterEfficiency=cms.untracked.double(1.0),
    pythiaHepMCVerbosity=cms.untracked.bool(False),
    comEnergy=cms.double(_args.comEnergy),
    PythiaParameters=cms.PSet(
        processParameters=cms.vstring(
            'Top:all = on',
        ),
        parameterSets=cms.vstring('processParameters'),
    ),
)

# Shared input tags, kept as cms.InputTag so the framework resolves them.
_genEventHepMC2 = cms.InputTag("generator", "unsmeared")
_genEventHepMC3 = cms.InputTag("generator", "unsmeared")
_genEventInfo = cms.InputTag("generator")

# --- (1) Standard path: GEN HepMC3 record via the reco::GenParticles ------
# The GenParticleProducer turns the HepMC2 record from the generator into a
# reco::GenParticle collection, exactly as a real GEN+SIM job would do. The
# GenParticles2HepMCConverter then rebuilds a HepMC3 GenEvent from those
# reco::GenParticles; this is the HepMC3 record Rivet has been consuming.
process.genParticles = cms.EDProducer("GenParticleProducer",
    abortOnUnknownPDGCode=cms.untracked.bool(False),
    saveBarCodes=cms.untracked.bool(True),
    src=cms.InputTag("generator:unsmeared"),
)
process.load("GeneratorInterface.RivetInterface.genParticles2HepMC_cfi")
process.genParticles2HepMC.genParticles = cms.InputTag("genParticles")
process.genParticles2HepMC.genEventInfo = _genEventInfo

# --- (2) Truth-graph path: GEN-only truth::Graph -> HepMC3 -----------------
# TruthGraphProducer requires SimTracks/SimVertices via evt.get(), which
# would throw if the products are absent. In a generator-only job there is no
# Geant4 step, so we emit empty SimTrack/SimVertex containers; the producer
# then builds a GEN-only raw TruthGraph (it handles an empty SimTrack/
# SimVertex container perfectly: the SIM-side loops are skipped).
process.emptySimContainers = cms.EDProducer("EmptySimContainersProducer",
    simTracks=cms.InputTag("emptySimContainers", "g4SimHits"),
    simVertices=cms.InputTag("emptySimContainers", "g4SimHits"),
)
_simTracksTag = cms.InputTag("emptySimContainers", "g4SimHits")
_simVerticesTag = cms.InputTag("emptySimContainers", "g4SimHits")

process.truthGraphProducer = cms.EDProducer(
    "TruthGraphProducer",
    genEventHepMC3=_genEventHepMC3,
    genEventHepMC=_genEventHepMC2,
    simTracks=_simTracksTag,
    simVertices=_simVerticesTag,
    addGenToSimEdges=cms.bool(True),   # no SimTracks => no GenToSim edges are produced
    collapseGenShower=cms.bool(False), # keep every GEN particle so the truth-graph
                                       # HepMC3 record matches the standard one as
                                       # closely as possible
)

# TruthLogicalGraphProducer reads the raw TruthGraph + (optional) SimTracks/
# SimVertices + HepMC2/HepMC3 to build the user-facing logical truth::Graph.
process.truthLogicalGraphProducer = cms.EDProducer(
    "TruthLogicalGraphProducer",
    src=cms.InputTag("truthGraphProducer"),
    simTracks=_simTracksTag,
    simVertices=_simVerticesTag,
    genEventHepMC3=_genEventHepMC3,
    genEventHepMC=_genEventHepMC2,
    rawGenPayload=cms.InputTag(""),  # only used during mixing, none here
    mergeGenSimVertices=cms.bool(True),
    postProcessing=cms.PSet(
        # Disable the post-processor's intermediate-resonance collapse so the
        # truth-graph record reproduces every GEN particle of the record - i.e.
        # the same set the standard GenParticles2HepMCConverter emits. The
        # expected differences then reduce to (a) the real beam protons (kept
        # by the truth graph, synthesized as dummies by the standard converter)
        # and (b) the vertex time (preserved by the truth graph, lost by the
        # reco::GenParticle), both of which leave the Rivet plots untouched.
        collapseIntermediateGenParticles=cms.bool(False),
        dropHitlessSimSubgraphs=cms.bool(False),  # GEN-only: no rechits
        seedPdgIds=cms.vint32(),
        seedHadronFlavors=cms.vint32(),
        seedParentDepth=cms.uint32(0),
        keepStableSpectators=cms.bool(True),
        attachSelectionSources=cms.bool(True),
        keepProductionSiblings=cms.bool(False),
        signalOnly=cms.bool(False),
        keepBunchCrossings=cms.vint32(),
        decayPdgIdGroups=cms.VPSet(),
        ignoredPdgIds=cms.vint32(),
        ignoredParticleIds=cms.vuint32(),
    ),
)

process.load("GeneratorInterface.RivetInterface.truthGraph2HepMC_cfi")
process.truthGraph2HepMC.truthGraph = cms.InputTag("truthLogicalGraphProducer")
process.truthGraph2HepMC.genEventInfo = _genEventInfo

# --- Two Rivet analyzers, one per path, with identical Rivet settings ----
process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

# Standard-path Rivet analyzer.
process.rivetAnalyzerStd = process.rivetAnalyzer.clone(
    AnalysisNames=cms.vstring(*_analysisNames),
    HepMCCollection=cms.InputTag("genParticles2HepMC:unsmeared"),
    GenEventInfoCollection=_genEventInfo,
    genLumiInfo=_genEventInfo,
    OutputFile=cms.string(_args.std_out),
    # Match the defaults of the cfi so the two paths are configured identically.
    useLHEweights=cms.bool(False),
    weightCap=cms.double(0.),
    NLOSmearing=cms.double(0.),
    skipMultiWeights=cms.bool(False),
    setIgnoreBeams=cms.bool(False),
    selectMultiWeights=cms.string(''),
    deselectMultiWeights=cms.string(''),
    setNominalWeightName=cms.string(''),
    CrossSection=cms.double(-1),
    DoFinalize=cms.bool(True),
)

# Truth-graph-path Rivet analyzer. Same Rivet settings; only the HepMC3
# input source differs (the truth-graph converter instead of the standard
# converter).
process.rivetAnalyzerTruthGraph = process.rivetAnalyzerStd.clone(
    HepMCCollection=cms.InputTag("truthGraph2HepMC:unsmeared"),
    OutputFile=cms.string(_args.truth_out),
)

# --- Path ----------------------------------------------------------------
# Run the generator first, then both parallel conversion chains, then the two
# Rivet analyzers. The two chains are independent: the standard chain is
# genParticles * genParticles2HepMC; the truth-graph chain is
# emptySimContainers * truthGraphProducer * truthLogicalGraphProducer *
# truthGraph2HepMC. The two Rivet analyzers follow. Since they share the same
# generator and the same Rivet analysis configuration, their YODA outputs are
# directly comparable.
process.p = cms.Path(
    process.generator
    * process.genParticles
    * process.genParticles2HepMC
    * process.emptySimContainers
    * process.truthGraphProducer
    * process.truthLogicalGraphProducer
    * process.truthGraph2HepMC
    * process.rivetAnalyzerStd
    * process.rivetAnalyzerTruthGraph
)

# Quiet the message logger, but keep the RivetAnalyzer and the new converter
# visible so a failure in either path is reported.
process.MessageLogger.cerr.threshold = cms.untracked.string("WARNING")
process.MessageLogger.cerr.RivetAnalyzer = cms.untracked.PSet(limit=cms.untracked.int32(-1))
process.MessageLogger.cerr.TruthGraph2HepMCConverter = cms.untracked.PSet(limit=cms.untracked.int32(-1))
process.MessageLogger.cerr.TruthLogicalGraphProducer = cms.untracked.PSet(limit=cms.untracked.int32(0))
process.MessageLogger.cerr.TruthGraphProducer = cms.untracked.PSet(limit=cms.untracked.int32(0))

print("=" * 78)
print("Rivet equivalence test (standard vs truth-graph converter)")
print("  Pythia8 ttbar at %g GeV, %d events" % (_args.comEnergy, _args.maxEvents))
print("  Rivet analyses: %s" % ", ".join(_analysisNames))
print("  Standard  path YODA: %s" % _args.std_out)
print("  Truth     path YODA: %s" % _args.truth_out)
print("  Compare with: python3 test/compareYodaFiles.py %s %s" % (_args.std_out, _args.truth_out))
print("=" * 78)
