# Pythia8DeadConeHooks_cfi.py
#
# Example configuration showing how to enable the b-quark dead-cone
# modification user hooks (MassRestorationHook or MasslessBWeightHook)
# in a CMS Pythia8 hadronizer job.
#
# The hooks are registered with the CustomHookFactory of the
# GeneratorInterface/Pythia8Interface package and are attached to Pythia via
# the "UserCustomization" parameter of the Pythia8HadronizerFilter /
# Pythia8HepMC3HadronizerFilter EDFilter.
#
# Two methods are provided:
#   * MassRestorationHook : switch off the dead cone by running the shower
#                           with a massless b ("5:m0 = 0.0") and restoring the b
#                           mass just before hadronization.
#   * MasslessBWeightHook  : soften the dead cone by reweighting trial
#                           emissions inside the dead cone. The b is kept massive
#                           in the shower. Requires the enhancement machinery to be
#                           switched on:
#                               Enhancements:doEnhanceTrial = on
#                               Enhancements:overSampleFSR  = 10
#
# The b mass used by the hooks defaults to 4.8 GeV and can be overridden via the
# "bMass" entry of the PSet.
#
# Usage: import this module (or copy the relevant blocks) into your generator
# configuration and select which hook to attach by adjusting the
# "UserCustomization" VPSet. Only one of the two methods should be enabled at a
# time.

import FWCore.ParameterSet.Config as cms

# ---------------------------------------------------------------------------
# Pythia8 settings blocks for the two dead-cone methods.
# ---------------------------------------------------------------------------

# Settings common to both methods (use the CP5 tune as an example; replace
# with the tune appropriate for your sample).
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

# Method 1: Mass Restoration.
# The b is showered massless ("5:m0 = 0.0") and the MassRestorationHook
# restores the b mass before hadronization.
pythia8DeadConeMassRestorationSettings = cms.vstring(
    '5:m0 = 0.0',
)

# Method 2: Massless b reweighting.
# The b keeps its mass in the shower; the trial emission probability is
# reweighted to enhance emissions inside the dead cone. The enhancement
# machinery must be switched on so the over-sampling actually produces the
# enhanced emissions.
pythia8DeadConeMasslessBWeightSettings = cms.vstring(
    'Enhancements:doEnhanceTrial = on',
    'Enhancements:overSampleFSR  = 10',
)

# ---------------------------------------------------------------------------
# Example generator block (Pythia8HadronizerFilter).
#
# This is an example for an internal-Pythia8 ttbar sample. Adapt the
# process, beam, energy, and tune to your needs. Pick the
# "UserCustomization" entry that corresponds to the method you want to use.
# ---------------------------------------------------------------------------
exampleDeadConeGenerator = cms.EDFilter("Pythia8HadronizerFilter",
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13000.),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        processParameters = cms.vstring(
            'Top:all = on',
        ),
        parameterSets = cms.vstring(
            'pythia8CommonSettings',
            'pythia8CP5Settings',
            'processParameters',
        )
    ),

    # -------------------------------------------------------------------
    # Pick ONE of the two UserCustomization blocks below.
    # -------------------------------------------------------------------

    # Method 1: Mass Restoration. Remember to also add the
    # pythia8DeadConeMassRestorationSettings block (containing '5:m0 = 0.0')
    # to the parameterSets above so the shower runs with a massless b.
    UserCustomization = cms.VPSet(
        cms.PSet(
            pluginName = cms.string("MassRestorationHook"),
            bMass = cms.double(4.8)
        )
    ),

    # Method 2: Massless b reweighting. Remember to also add the
    # pythia8DeadConeMasslessBWeightSettings block (containing the
    # Enhancements:* settings) to the parameterSets above.
    # UserCustomization = cms.VPSet(
    #     cms.PSet(
    #         pluginName = cms.string("MasslessBWeightHook"),
    #         bMass = cms.double(4.8)
    #     )
    # ),
)
