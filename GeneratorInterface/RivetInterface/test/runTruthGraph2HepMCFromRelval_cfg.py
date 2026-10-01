# Minimal end-to-end run of the new TruthGraph2HepMCConverter on a RECO relval
# file (step3.root) that already carries both graph products from the in-time
# truth graph step:
#   - "mix":TruthGraph                    "HLT"  (raw TruthGraph)
#   - "truthLogicalGraphProducer":truth::Graph "HLT"  (logical truth::Graph)
# So this config does not need to run the truth-graph producers: it only reads
# the existing truth::Graph and converts it to an edm::HepMC3Product, the way a
# Rivet/PARTICLE-LEVEL workflow will. The HepMC3 record is also written to a
# text file for inspection.
#
# Usage:
#   cmsRun test/runTruthGraph2HepMCFromRelval_cfg.py [INPUT_ROOT] [-n N]
# Defaults to the relval step3.root (RECO) produced with runTheMatrix.

import os
import argparse
import FWCore.ParameterSet.Config as cms

_parser = argparse.ArgumentParser(description="TruthGraph2HepMCConverter on relval step3")
_parser.add_argument("inputRoot", nargs="?",
                    default=os.path.join(
                        os.environ.get("CMSSW_BASE"),
                        "src/38434.0_TTbar_14TeV+Run4D128/step3.root"),
                    help="input relval file (default: %(default)s)")
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
_out_hepmc = os.path.join(_outdir, "truthgraph2HepMC.fromRelval.hepmc")

process = cms.Process("TRUTH2HEPMC")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(_args.maxEvents))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring("file:" + _input_path))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

# The ParticleDataTable is consumed by the converter to look up the generated
# mass by PDG id, matching what the standard converter does.
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load("GeneratorInterface.RivetInterface.truthGraph2HepMC_cfi")
# The truth::Graph was produced by the in-time truth-graph step on the file:
process.truthGraph2HepMC.truthGraph = cms.InputTag("truthLogicalGraphProducer", "", "HLT")
# The generator's run/event info (signalProcessID, qScale, alphaQCD/QED, weights,
# PDF, cross-section) still comes from the generatorSmeared module label.
process.truthGraph2HepMC.genEventInfo = cms.InputTag("generator")
process.truthGraph2HepMC.writeHepMC = cms.untracked.bool(True)
process.truthGraph2HepMC.outputFile = cms.untracked.string(_out_hepmc)
process.truthGraph2HepMC.verbosity = cms.untracked.uint32(1)

process.p = cms.Path(process.truthGraph2HepMC)

process.MessageLogger.cerr.threshold = cms.untracked.string("WARNING")
process.MessageLogger.cerr.TruthGraph2HepMCConverter = cms.untracked.PSet(limit=cms.untracked.int32(-1))

print("=" * 78)
print("TruthGraph2HepMCConverter on relval step3 (uses already-produced truth graph)")
print("  input       :", _input_path)
print("  events      :", _args.maxEvents)
print("  HepMC3 file :", _out_hepmc)
print("=" * 78)
