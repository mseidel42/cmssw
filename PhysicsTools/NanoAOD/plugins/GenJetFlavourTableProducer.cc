// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
class GenJetFlavourTableProducer : public edm::stream::EDProducer<> {
public:
  explicit GenJetFlavourTableProducer(const edm::ParameterSet& iConfig)
      : name_(iConfig.getParameter<std::string>("name")),
        src_(consumes<std::vector<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("src"))),
        cut_(iConfig.getParameter<std::string>("cut"), true),
        deltaR_(iConfig.getParameter<double>("deltaR")),
        jetFlavourInfosToken_(
            consumes<reco::JetFlavourInfoMatchingCollection>(iConfig.getParameter<edm::InputTag>("jetFlavourInfos"))) {
    produces<nanoaod::FlatTable>();
  }

  ~GenJetFlavourTableProducer() override {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src")->setComment("input genJet collection");
    desc.add<edm::InputTag>("jetFlavourInfos")->setComment("input flavour info collection");
    desc.add<std::string>("name")->setComment("name of the genJet FlatTable we are extending with flavour information");
    desc.add<std::string>("cut")->setComment("cut on input genJet collection");
    desc.add<double>("deltaR")->setComment("deltaR to match genjets");
    descriptions.add("genJetFlavourTable", desc);
  }

private:
  void produce(edm::Event&, edm::EventSetup const&) override;

  uint32_t fjContribFlavArrayToInt(const std::vector<int>& fjContribFlav) const;
  int16_t  fjContribFlavArrayToLeading(const std::vector<int>& fjContribFlav) const;

  std::string name_;
  edm::EDGetTokenT<std::vector<reco::GenJet> > src_;
  const StringCutObjectSelector<reco::GenJet> cut_;
  const double deltaR_;
  edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavourInfosToken_;
};

// ------------ method called to produce the data  ------------
void GenJetFlavourTableProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const auto& jetsProd = iEvent.get(src_);
  const auto& jetFlavourInfosProd = iEvent.get(jetFlavourInfosToken_);

  unsigned int ncand = 0;
  std::vector<int16_t> partonFlavour;
  std::vector<uint8_t> hadronFlavour;
  std::vector<uint8_t> nBHadrons;
  std::vector<uint8_t> nCHadrons;
  // fastjet::contrib flavour info
  std::vector<uint32_t> fjContribGHSAlgoFlav;
  std::vector<uint32_t> fjContribGHSFullAlgoFlav;
  std::vector<int16_t> fjContribGHSAlgoLeadingFlav;
  std::vector<int16_t> fjContribGHSFullAlgoleadingFlav;


  for (const reco::GenJet& jet : jetsProd) {
    if (!cut_(jet))
      continue;
    ++ncand;
    bool matched = false;
    for (const reco::JetFlavourInfoMatching& jetFlavourInfoMatching : jetFlavourInfosProd) {
      if (deltaR(jet.p4(), jetFlavourInfoMatching.first->p4()) < deltaR_) {
        partonFlavour.push_back(jetFlavourInfoMatching.second.getPartonFlavour());
        hadronFlavour.push_back(jetFlavourInfoMatching.second.getHadronFlavour());
        nBHadrons.push_back(jetFlavourInfoMatching.second.getbHadrons().size());
        nCHadrons.push_back(jetFlavourInfoMatching.second.getcHadrons().size());
        // fastjet::contrib flavour info
        if (jetFlavourInfoMatching.second.haveAlgoFlav(reco::FlavAlgo::kGHS)) {
          fjContribGHSAlgoFlav.push_back(
            fjContribFlavArrayToInt(jetFlavourInfoMatching.second.getAlgoFlav(reco::FlavAlgo::kGHS)));
          fjContribGHSAlgoLeadingFlav.push_back(
            fjContribFlavArrayToLeading(jetFlavourInfoMatching.second.getAlgoFlav(reco::FlavAlgo::kGHS)));
        } else {
          fjContribGHSAlgoFlav.push_back(0);
          fjContribGHSAlgoLeadingFlav.push_back(0);
        }
        if (jetFlavourInfoMatching.second.haveAlgoFlav(reco::FlavAlgo::kGHSFull)) {
          fjContribGHSFullAlgoFlav.push_back(
            fjContribFlavArrayToInt(jetFlavourInfoMatching.second.getAlgoFlav(reco::FlavAlgo::kGHSFull)));
          fjContribGHSFullAlgoleadingFlav.push_back(
            fjContribFlavArrayToLeading(jetFlavourInfoMatching.second.getAlgoFlav(reco::FlavAlgo::kGHSFull)));
        } else {
          fjContribGHSFullAlgoFlav.push_back(0);
          fjContribGHSFullAlgoleadingFlav.push_back(0);
        }
        matched = true;
        break;
      }
    }
    if (!matched) {
      partonFlavour.push_back(0);
      hadronFlavour.push_back(0);
      nBHadrons.push_back(0);
      nCHadrons.push_back(0);
      // fastjet::contrib flavour info
      fjContribGHSAlgoFlav.push_back(0);
      fjContribGHSFullAlgoFlav.push_back(0);
      fjContribGHSAlgoLeadingFlav.push_back(0);
      fjContribGHSFullAlgoleadingFlav.push_back(0);
    }
  }

  auto tab = std::make_unique<nanoaod::FlatTable>(ncand, name_, false, true);
  tab->addColumn<int16_t>("partonFlavour", partonFlavour, "flavour from parton matching");
  tab->addColumn<uint8_t>("hadronFlavour", hadronFlavour, "flavour from hadron ghost clustering");
  tab->addColumn<uint8_t>("nBHadrons", nBHadrons, "number of b-hadrons");
  tab->addColumn<uint8_t>("nCHadrons", nCHadrons, "number of c-hadrons");
  tab->addColumn<uint32_t>("AlgoGHSFlav", fjContribGHSAlgoFlav,
                           "fastjet::contrib GHS flavour info, encoded as integer");
  tab->addColumn<uint32_t>("AlgoGHSFullFlav", fjContribGHSFullAlgoFlav,
                            "fastjet::contrib GHS full algorithm flavour info, encoded as integer");
  tab->addColumn<int16_t>("AlgoGHSLeadFlav", fjContribGHSAlgoLeadingFlav,
                            "fastjet::contrib GHS flavour info, leading flavour only, encoded as integer");
  tab->addColumn<int16_t>("AlgoGHSFullLeadFlav", fjContribGHSFullAlgoleadingFlav,
                            "fastjet::contrib GHS full algorithm flavour info, leading flavour only, encoded as integer");               

  iEvent.put(std::move(tab));
}

/// ------------- Convert fastjet::contrib flavour algo result from array to integer -------------
uint32_t GenJetFlavourTableProducer::fjContribFlavArrayToInt(const std::vector<int>& fjContribFlav) const {
  /// Coding rule: lowest 3 bits for 0th element from the array, the flag;
  ///              then 4 bits from lowest to highest, in the order "d, u, s, c, b, t"
  ///              For each flavour, highest bit 0/1 is q/qbar,
  ///              and next 3 bits are "none, 1, 2, 3, 4, 5, even and >=6, odd and >=7"
  uint32_t result = static_cast<uint32_t>(fjContribFlav[0] & 0x7);
  uint32_t flavAbs = 0;
  for (unsigned int i = 1; i < fjContribFlav.size(); ++i) {
    flavAbs = static_cast<uint32_t>(fjContribFlav[i] > 0 ? fjContribFlav[i] : -fjContribFlav[i]);
    if (flavAbs > 7)
      flavAbs = (6 | ((fjContribFlav[i] & 1))); // even or odd and >= 6
    // Handle the sign.
    if (fjContribFlav[i] < 0)
      flavAbs |= (1 << 3); // set the sign bit for qbar
    result |= (flavAbs << (4 * i - 1));
  }
  return result;
}

int16_t GenJetFlavourTableProducer::fjContribFlavArrayToLeading(const std::vector<int>& fjContribFlav) const {
  /// Taking only the heaviest flavour from the array. q/qbar separated by +/-.
  /// Mark only the existence of the flavour, not the number of constituents.
  int16_t result = 0;
  for (unsigned int i = fjContribFlav.size() - 1; i > 0; --i) {
    if (fjContribFlav[i] != 0) {
      result = static_cast<int16_t>(fjContribFlav[i] > 0 ? i : -i);
      break;
    }
  }
  return result;
}

#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE(GenJetFlavourTableProducer);
