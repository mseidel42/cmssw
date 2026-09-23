#ifndef GeneratorInterface_Pythia8Interface_MassRestoration_h
#define GeneratorInterface_Pythia8Interface_MassRestoration_h

// MassRestoration.h
//
// CMSSW integration of the MassRestorationHook developed for the standalone
// Pythia8 main115.cc example that switches off the b-quark dead cone.
//
// The hook is to be used together with a massless b in the shower, i.e. the
// Pythia setting
//     5:m0 = 0.0
// must be set in the configuration. After the shower, but before hadronization,
// the b mass is restored. For each final-state b quark the nearest (smallest
// opening angle) gluon is found such that the pair invariant mass
// s > (m_b + m_g)^2. The b and gluon are then put on their new mass shells,
// keeping the pair 3-momentum fixed in the pair rest frame. If no
// kinematically allowed partner can be found for a b quark the event is
// vetoed (Pythia re-tries the parton level).
//
// The b mass used for the restoration is configurable via the ParameterSet
// entry "bMass" (defaults to 4.8 GeV, the typical B-hadron mass).

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "Pythia8/Pythia.h"
#include <memory>
#include "GeneratorInterface/Pythia8Interface/interface/CustomHook.h"

class MassRestorationHook : public Pythia8::UserHooks {
public:
  MassRestorationHook(const edm::ParameterSet& iConfig);
  ~MassRestorationHook() override {}

  // Veto at the end of the parton level, just before hadronization.
  bool canVetoPartonLevel() override { return true; }
  // Modify the event record (restore b mass) and let hadronization proceed.
  bool doVetoPartonLevel(const Pythia8::Event& constEvent) override;

private:
  double bMass_;
};

REGISTER_USERHOOK(MassRestorationHook);
#endif
