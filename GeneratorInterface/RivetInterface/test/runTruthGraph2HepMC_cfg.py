# Standalone end-to-end run of the new TruthGraph2HepMCConverter on the relval
# GEN+SIM file. Produces the raw TruthGraph from the HepMC2 record and the
# SimTrack/SimVertex collections, builds the logical truth::Graph, and converts
# it to a HepMC3 GenEvent that can be fed to Rivet (RivetAnalyzer /
# ParticleLevelProducer). The HepMC3 record is also written to a text file
# (events.hepmc) so the user can inspect the output by eye or diff it against the
# standard converter's output.
#
# Usage:
#   cmsRun test/runTruthGraph2HepMC_cfg.py [INPUT_ROOT] [-n N]
# Defaults to the relval step1.root (GEN+SIM) produced with runTheMatrix.
#
# For a self-contained equivalence test against the standard converter, see
# compareTruthGraphVsGenParticles_cfg.py instead.

import os
import argparse
import FWCore.ParameterSet.Config as cms

_parser = argparse.ArgumentParser(description="TruthGraph2HepMCConverter standalone run")
_parser.add_argument("inputRoot", nargs="?",
                    default=os.path.join(
                        os.environ.get("CMSSW_BASE"),
                        "src/38434.0_TTbar_14TeV+Run4D128/step1.root"),
                    help="input GEN+SIM relval file (default: %(default)s)")
_parser.add_argument("-n", "--maxEvents", type=int, default=2,
                    help="number of events to process (default: %(default)s)")
_args, _unknown = _parser.parse_known_args()

_input_path = _args.inputRoot
if not os.path.isabs(_input_path) and not os.path.exists(_input_path):
    _here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
    _try = os.path.join(_here, _input_path)
    if os.path.exists(_try):
        _input_path = _try
if not os.path.exists(_input_path):
    raise RuntimeError("Input file not found: %s" % _input_path)

_outdir = os.path.dirname(_input_path) or "."
_out_hepmc = os.path.join(_outdir, "truthgraph2HepMC.events.hepmc")

process = cms.Process("TRUTH2HEPMC")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(_args.maxEvents))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring("file:" + _input_path))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

# The ParticleDataTable is consumed by the converter to look up the generated
# mass by PDG id, matching what the standard converter does.
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# TruthGraphProducer: builds the raw TruthGraph (the GEN+SIM merged record) from
# the HepMC2 generator record (HepMC3 may also be present in the file) plus the
# SimTrack/SimVertex collections.
process.truthGraphProducer = cms.EDProducer(
    "TruthGraphProducer",
    genEventHepMC3=cms.InputTag("generatorSmeared"),
    genEventHepMC=cms.InputTag("generatorSmeared"),
    simTracks=cms.InputTag("g4SimHits"),
    simVertices=cms.InputTag("g4SimHits"),
    addGenToSimEdges=cms.bool(True),
    collapseGenShower=cms.bool(False),
)

# TruthLogicalGraphProducer: builds the user-facing logical truth::Graph from the
# raw graph (this is what the converter reads). The post-processor is set to the
# defaults that keep the full graph: no selection, no pile-up filtering, so the
# truth-graph HepMC3 record carries the entire GEN record exactly as the
# standard converter sees it.
process.truthLogicalGraphProducer = cms.EDProducer(
    "TruthLogicalGraphProducer",
    src=cms.InputTag("truthGraphProducer"),
    simTracks=cms.InputTag("g4SimHits"),
    simVertices=cms.InputTag("g4SimHits"),
    genEventHepMC3=cms.InputTag("generatorSmeared"),
    genEventHepMC=cms.InputTag("generatorSmeared"),
    rawGenPayload=cms.InputTag(""),
    mergeGenSimVertices=cms.bool(True),
    postProcessing=cms.PSet(
        collapseIntermediateGenParticles=cms.bool(True),
        dropHitlessSimSubgraphs=cms.bool(False),
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
process.truthGraph2HepMC.genEventInfo = cms.InputTag("generator")
process.truthGraph2HepMC.writeHepMC = cms.untracked.bool(True)
process.truthGraph2HepMC.outputFile = cms.untracked.string(_out_hepmc)
process.truthGraph2HepMC.verbosity = cms.untracked.uint32(1)

process.p = cms.Path(
    process.truthGraphProducer
    * process.truthLogicalGraphProducer
    * process.truthGraph2HepMC
)

process.MessageLogger.cerr.threshold = cms.untracked.string("WARNING")
process.MessageLogger.cerr.TruthGraph2HepMCConverter = cms.untracked.PSet(limit=cms.untracked.int32(-1))
process.MessageLogger.cerr.TruthLogicalGraphProducer = cms.untracked.PSet(limit=cms.untracked.int32(0))
process.MessageLogger.cerr.TruthGraphProducer = cms.untracked.PSet(limit=cms.untracked.int32(0))

print("=" * 78)
print("TruthGraph2HepMCConverter standalone run")
print("  input       :", _input_path)
print("  events     :", _args.maxEvents)
print("  HepMC3 file:", _out_hepmc)
print("=" * 78)
