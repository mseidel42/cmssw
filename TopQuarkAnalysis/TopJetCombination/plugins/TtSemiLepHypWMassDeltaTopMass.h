#ifndef TtSemiLepHypWMassDeltaTopMass_h
#define TtSemiLepHypWMassDeltaTopMass_h

#include "TopQuarkAnalysis/TopJetCombination/interface/TtSemiLepHypothesis.h"

#include "DataFormats/PatCandidates/interface/Particle.h"

class TtSemiLepHypWMassDeltaTopMass : public TtSemiLepHypothesis  {

 public:

  explicit TtSemiLepHypWMassDeltaTopMass(const edm::ParameterSet& cfg);
  ~TtSemiLepHypWMassDeltaTopMass() override {};

 private:

  /// build the event hypothesis key
  void buildKey() override { key_= TtSemiLeptonicEvent::kWMassDeltaTopMass; };  
  /// build event hypothesis from the reco objects of a semi-leptonic event 
  void buildHypo(edm::Event& evt,
			 const edm::Handle<edm::View<reco::RecoCandidate> >& leps,
			 const edm::Handle<std::vector<pat::MET> >& mets,
			 const edm::Handle<std::vector<pat::Jet> >& jets,
			 std::vector<int>& match, const unsigned int iComb) override; // { TtSemiLepHypothesis::buildHypo(leps, mets, jets, match); };

//   edm::EDGetTokenT<std::vector<int> > statusToken_;
//   edm::EDGetTokenT<std::vector<pat::Particle> > partonsHadPToken_;
//   edm::EDGetTokenT<std::vector<pat::Particle> > partonsHadQToken_;
//   edm::EDGetTokenT<std::vector<pat::Particle> > partonsHadBToken_;
//   edm::EDGetTokenT<std::vector<pat::Particle> > partonsLepBToken_;
//   edm::EDGetTokenT<std::vector<pat::Particle> > leptonsToken_;
  edm::EDGetTokenT<std::vector<pat::Particle> > neutrinosToken_;

};

#endif
