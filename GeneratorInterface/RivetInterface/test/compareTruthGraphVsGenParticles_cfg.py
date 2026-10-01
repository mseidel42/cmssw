# End-to-end equivalence test for the new TruthGraph2HepMCConverter against the
# standard GenParticles2HepMCConverter. The new "Truth Graph" feature produces a
# truth::Graph (from GEN+SIM, via TruthLogicalGraphProducer); this config feeds
# the same record to the truth-graph converter, and runs the standard GEN-only
# converter from the same reco::GenParticles, on the SAME event, then writes
# both HepMC3 text files side-by-side for a physics-equivalence diff.
#
# The two HepMC3 records produced here are NOT expected to be byte-identical:
#   - The truth-graph record carries the two real beam protons from the HepMC
#     generator record (status 4), reproducing the actual hard-scatter vertex; the
#     standard converter can only synthesize dummy protons since the reco::GenParticles
#     collection dropped the beam particles.
#   - The truth-graph record preserves the generator's vertex time (the graph
#     stores ns, converted to mm-of-c*t); the standard converter sets t=0 since
#     the reco::GenParticle does not carry a vertex time.
# The physics content of the GEN record (the particles Rivet analyses consume) is
# identical between the two. Use the compareHepMCFiles.py helper to verify the
# equivalence instead of a plain `diff`.
#
# Usage:
#   cmsRun test/compareTruthGraphVsGenParticles_cfg.py [INPUT_ROOT] [-n N]
# Defaults to the relval step1.root (GEN+SIM) produced with runTheMatrix.
# Run from the test/ directory or from the src/ area; the input path resolves
# relative to it.

import os
import sys
import argparse
import FWCore.ParameterSet.Config as cms

# ---- Argument parsing (done first, in cmsRun-friendly way) ----------------
_parser = argparse.ArgumentParser(description="TruthGraph2HepMC equivalence test")
# Take one positional INPUT_ROOT plus optional -n.
_parser.add_argument("inputRoot", nargs="?",
                    default=os.path.join(
                        os.environ.get("CMSSW_BASE"),
                        "src/38434.0_TTbar_14TeV+Run4D128/step1.root"),
                    help="input GEN+SIM relval file (default: %(default)s)")
_parser.add_argument("-n", "--maxEvents", type=int, default=2,
                    help="number of events to process (default: %(default)s)")
_args, _unknown = _parser.parse_known_args()

# Resolve the input relative to the test/ directory if it looks relative and the
# raw path does not exist. This makes the test work when cmsRun is run from
# src/GeneratorInterface/RivetInterface/test/.
_input_path = _args.inputRoot
if not os.path.isabs(_input_path) and not os.path.exists(_input_path):
    _here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
    _try = os.path.join(_here, _input_path)
    if os.path.exists(_try):
        _input_path = _try

if not os.path.exists(_input_path):
    raise RuntimeError("Input file not found: %s" % _input_path)

_outdir = os.path.dirname(_input_path) or "."
_out_truth = os.path.join(_outdir, "truthgraph2HepMC.events.hepmc")
_out_std = os.path.join(_outdir, "genParticles2HepMC.events.hepmc")

# ---- The process ---------------------------------------------------------
process = cms.Process("TRUTH2HEPMC")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(_args.maxEvents))

process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring("file:" + _input_path))

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

# The particle data table is needed by both converters to look up the generated
# mass by PDG id (so both HepMC3 records carry the same mass for a given particle).
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# Standard track/sim-vertex collections: needed by the TruthGraphProducer so the
# raw truth graph merges in the SIM side. On the GEN+SIM relval step they live
# under the SIM process and the g4SimHits module label.
_truthGraphGenHepMC2 = cms.InputTag("generatorSmeared")
_truthGraphGenHepMC3 = cms.InputTag("generatorSmeared")
_simTracks = cms.InputTag("g4SimHits")
_simVertices = cms.InputTag("g4SimHits")

# --- Standard converter (the one Rivet was already using in CMSSW) ---------
# It reads the reco::GenParticle collection ("genParticles") and the
# GenEventInfoProduct ("generator"), producing an edm::HepMC3Product:unsmeared.
process.load("GeneratorInterface.RivetInterface.genParticles2HepMC_cfi")
# On the GEN+SIM relval step the genParticles live under the SIM process label.
process.genParticles2HepMC.genParticles = cms.InputTag("genParticles", "", "SIM")
process.genParticles2HepMC.genEventInfo = cms.InputTag("generator")
process.genParticles2HepMC.writeHepMC = cms.untracked.bool(True)
process.genParticles2HepMC.outputFile = cms.untracked.string(_out_std)

# --- TruthGraph producers (the new "Truth Graph" feature) ------------------
process.truthGraphProducer = cms.EDProducer(
    "TruthGraphProducer",
    genEventHepMC3=_truthGraphGenHepMC3,
    genEventHepMC=_truthGraphGenHepMC2,
    simTracks=_simTracks,
    simVertices=_simVertices,
    addGenToSimEdges=cms.bool(True),
    collapseGenShower=cms.bool(False),  # keep the full GEN record: matches the
                                        # reco::GenParticles the standard converter reads
)

# The post-processor default keeps the full graph (seedPdgIds = {}), so the
# truth-graph record carries every GEN particle of the record (signal + pileup)
# exactly as the standard converter sees them.
process.truthLogicalGraphProducer = cms.EDProducer(
    "TruthLogicalGraphProducer",
    src=cms.InputTag("truthGraphProducer"),
    simTracks=_simTracks,
    simVertices=_simVertices,
    genEventHepMC3=_truthGraphGenHepMC3,
    genEventHepMC=_truthGraphGenHepMC2,
    rawGenPayload=cms.InputTag(""),  # only used in mixing, none here
    mergeGenSimVertices=cms.bool(True),
    postProcessing=cms.PSet(
        # Disable the post-processor's "collapse intermediate GEN particles" rule
        # so the truth-graph HepMC3 record reproduces every GEN particle of the
        # record - i.e. the same set the standard GenParticles2HepMCConverter
        # emits. With this off, the only expected differences vs. the standard
        # converter are the real beam protons (kept here) vs. the standard
        # converter's dummy beam protons, and the vertex time, which the
        # truth graph preserves and the reco::GenParticle does not carry.
        collapseIntermediateGenParticles=cms.bool(False),
        dropHitlessSimSubgraphs=cms.bool(False),  # GEN+SIM relval step1: no rechits yet
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

# --- New converter: truth::Graph -> HepMC3 ---------------------------------
process.load("GeneratorInterface.RivetInterface.truthGraph2HepMC_cfi")
process.truthGraph2HepMC.truthGraph = cms.InputTag("truthLogicalGraphProducer")
process.truthGraph2HepMC.genEventInfo = cms.InputTag("generator")
process.truthGraph2HepMC.writeHepMC = cms.untracked.bool(True)
process.truthGraph2HepMC.outputFile = cms.untracked.string(_out_truth)
process.truthGraph2HepMC.verbosity = cms.untracked.uint32(1)

# --- Path ----------------------------------------------------------------
process.p = cms.Path(
    process.genParticles2HepMC
    * process.truthGraphProducer
    * process.truthLogicalGraphProducer
    * process.truthGraph2HepMC
)

process.MessageLogger.cerr.threshold = cms.untracked.string("WARNING")
process.MessageLogger.cerr.TruthGraph2HepMCConverter = cms.untracked.PSet(limit=cms.untracked.int32(-1))
process.MessageLogger.cerr.TruthLogicalGraphProducer = cms.untracked.PSet(limit=cms.untracked.int32(0))
process.MessageLogger.cerr.TruthGraphProducer = cms.untracked.PSet(limit=cms.untracked.int32(0))

# ---- Print what we ran, so the user can re-run the comparison easily ----
print("=" * 78)
print("TruthGraph2HepMC equivalence test")
print("  input        :", _input_path)
print("  events        :", _args.maxEvents)
print("  std HepMC3    :", _out_std)
print("  truth HepMC3  :", _out_truth)
print("  compare with  : python3 test/compareHepMCFiles.py %s %s" % (_out_std, _out_truth))
print("=" * 78)
