#ifndef GeneratorInterface_Pythia8Interface_MasslessLHEInput_h
#define GeneratorInterface_Pythia8Interface_MasslessLHEInput_h

// MasslessLHEInput.h
//
// CMSSW hook that massless-izes the b quarks coming from the LHE-input event
// record before the shower runs, so that the FSR shower sees no dead cone for
// them.
//
// Why this is needed:
// The Pythia setting "5:m0 = 0.0" only sets the *default* b mass used when a
// new b quark is created (e.g. by g -> bb inside the shower). It does NOT
// affect b quarks that are already in the event record from the LHE, which
// carry the mass that Powheg (or another ME generator) wrote into the LHE.
// Those b quarks are then showered with their LHE mass, i.e. with an active
// dead cone. To genuinely switch off the dead cone when showering massive
// LHE-input b quarks, this hook additionally massless-izes the final-state
// b quarks at the process level (after the LHE is read, before the shower).
//
// Top-mass preservation:
// When the b quark comes from a resonance decay (e.g. t -> b W in a Powheg
// LHE file with full spin-correlated top decays), the hook rescales all the
// resonance daughters by a common factor alpha in the resonance rest frame.
// alpha is chosen so that (i) the b becomes massless (E_b' = alpha |p_b|) and
// (ii) the total energy of the daughters equals the resonance mass. This
// preserves the resonance (top) mass exactly. For the two-body t -> b W case
// the formula reduces to the standard two-body kinematics with m_b = 0; for
// multi-body decays (e.g. bb4l-style top -> b W g) the uniform-rescaling
// approach is a reasonable approximation that preserves the resonance mass
// and the directions of all decay products.
//
// For b quarks that are not from a resonance decay (e.g. b-initiated hard
// processes in the 5FS, or b from g -> bb in the ME), the hook simply
// preserves the 3-momentum, sets E = |p|, and sets the stored mass to 0.
//
// Incoming b quarks are NOT modified -- changing their 4-momentum would
// distort the hard ME PDF kinematics. (Incoming b's have negative status and
// are not "final", so they are automatically excluded by the isFinal()
// check.)

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "Pythia8/Pythia.h"
#include <memory>
#include "GeneratorInterface/Pythia8Interface/interface/CustomHook.h"

class MasslessLHEInputHook : public Pythia8::UserHooks {
public:
  MasslessLHEInputHook(const edm::ParameterSet& iConfig);
  ~MasslessLHEInputHook() override {}

  // Veto at the process level, after the LHE is read and before the shower.
  bool canVetoProcessLevel() override { return true; }
  bool doVetoProcessLevel(Pythia8::Event& event) override;
};

REGISTER_USERHOOK(MasslessLHEInputHook);
#endif
