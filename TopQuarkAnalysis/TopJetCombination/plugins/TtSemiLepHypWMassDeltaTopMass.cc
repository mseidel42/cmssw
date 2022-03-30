#include "TopQuarkAnalysis/TopJetCombination/plugins/TtSemiLepHypWMassDeltaTopMass.h"

TtSemiLepHypWMassDeltaTopMass::TtSemiLepHypWMassDeltaTopMass(const edm::ParameterSet& cfg)
    : TtSemiLepHypothesis(cfg),
      neutrinosToken_(consumes<std::vector<pat::Particle> >(cfg.getParameter<edm::InputTag>("neutrinos"))) {}

void TtSemiLepHypWMassDeltaTopMass::buildHypo(edm::Event& evt,
                                              const edm::Handle<edm::View<reco::RecoCandidate> >& leps,
                                              const edm::Handle<std::vector<pat::MET> >& mets,
                                              const edm::Handle<std::vector<pat::Jet> >& jets,
                                              std::vector<int>& match,
                                              const unsigned int iComb) {
  //   edm::Handle<std::vector<int> > status;
  //   evt.getByToken(statusToken_, status);
  //   if( (*status)[iComb] != 0 ){
  //     // create empty hypothesis if kinematic fit did not converge
  //     return;
  //   }
  edm::Handle<std::vector<pat::Particle> > neutrinos;
  evt.getByToken(neutrinosToken_, neutrinos);

  // -----------------------------------------------------
  // add neutrino
  // -----------------------------------------------------

  if (!neutrinos->empty()) {
    setCandidate(neutrinos, iComb, neutrino_);
  }
}