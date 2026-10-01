import FWCore.ParameterSet.Config as cms

# TruthGraph2HepMCConverter is the truth-graph twin of the standard
# GenParticles2HepMCConverter: it reads the user-facing logical truth::Graph
# produced by TruthLogicalGraphProducer (the new "Truth Graph" feature that
# combines the GEN and SIM truth records) and rebuilds a HepMC3 GenEvent from
# it so the same record can be fed to Rivet (RivetAnalyzer / ParticleLevelProducer),
# exactly as the standard converter lets it be fed from the GEN reco::GenParticles.
#
# By default only the GEN side of the graph is emitted (matching the scope of the
# standard converter); the SIM-only Geant4 secondaries can be added with
# includeSimOnlyParticles=True.

truthGraph2HepMC = cms.EDProducer("TruthGraph2HepMCConverter",
    # The truth::Graph product (typically the output of TruthLogicalGraphProducer).
    truthGraph = cms.InputTag("truthLogicalGraphProducer"),
    # Generator module label that provides the GenEventInfoProduct AND the
    # GenRunInfoProduct (same InputTag, looked up per-event and per-run). This is
    # the same convention the standard GenParticles2HepMCConverter uses, so the
    # truth-graph and the standard HepMC3 records carry the same event attributes
    # (signal_process_id, qScale, alphaQCD, alphaQED, weights, pdf info, cross
    # section).
    genEventInfo = cms.InputTag("generator"),
)
