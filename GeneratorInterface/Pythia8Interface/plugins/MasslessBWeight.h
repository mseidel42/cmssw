#ifndef GeneratorInterface_Pythia8Interface_MasslessBWeight_h
#define GeneratorInterface_Pythia8Interface_MasslessBWeight_h

// MasslessBWeight.h
//
// CMSSW integration of the MasslessBWeightHook developed for the standalone
// Pythia8 main115.cc example that softens the b-quark dead cone.
//
// Unlike MassRestorationHook, the b is kept massive in the shower, but the
// trial emission probability is reweighted so that emissions inside the
// dead cone are enhanced. To be used together with the Pythia settings
//     Enhancements:doEnhanceTrial = on
//     Enhancements:overSampleFSR  = 10
// so the over-sampling actually generates the enhanced emissions.
//
// The reweighting factor is
//     W = (1 + (theta0/theta)^2)^2,   theta0 = m_b / E
// capped at 1000 to protect against numerical spikes in the deep collinear
// region. The b mass used in the weight is configurable via the ParameterSet
// entry "bMass" (defaults to 4.8 GeV, the typical B-hadron mass).

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "Pythia8/Pythia.h"
#include <memory>
#include "GeneratorInterface/Pythia8Interface/interface/CustomHook.h"

class MasslessBWeightHook : public Pythia8::UserHooks {
public:
  MasslessBWeightHook(const edm::ParameterSet& iConfig);
  ~MasslessBWeightHook() override {}

  // Tell Pythia we want to enhance trial emissions.
  bool canEnhanceEmission() override { return true; }
  // The core reweighting logic. iRad, iRec, iEmt are indices in workEvent.
  double enhanceEmission(int beamKind, int iRad, int iRec, int iEmt) override;

private:
  double mB;
};

REGISTER_USERHOOK(MasslessBWeightHook);
#endif
