#ifndef PAT_SMEAREDJETPRODUCERT_H
#define PAT_SMEAREDJETPRODUCERT_H

/** \class SmearedJetProducerT
 *
 * Produce collection of "smeared" jets.
 *
 * The aim of this correction is to account for the difference in jet energy resolution
 * between Monte Carlo simulation and Data.
 *
 * \author Sébastien Brochet
 *
 */

#include "CommonTools/Utils/interface/PtComparator.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>

#include <memory>
#include <random>

using std::endl;
using std::cout;

namespace reco {
  class PFJet;
  class CaloJet;
}

namespace pat {
  class Jet;
  class GenJetMatcher {
  public:
    GenJetMatcher(const edm::ParameterSet& cfg, edm::ConsumesCollector&& collector)
        : m_genJetsTag(cfg.getParameter<edm::InputTag>("genJets")),
          m_genJetsToken(collector.consumes<reco::GenJetCollection>(m_genJetsTag)),
          m_use_own_match(m_genJetsTag.label() == ""),
          m_dR_max(cfg.getParameter<double>("dRMax")),
          m_dPt_max_factor(cfg.getParameter<double>("dPtMaxFactor")) {
      // Empty
    }

    static void fillDescriptions(edm::ParameterSetDescription& desc) {
      desc.add<edm::InputTag>("genJets", edm::InputTag(""))->setComment("Empty: attempt using pat::Jet internal genJets.");
      desc.add<double>("dRMax", 0.2)->setComment("= cone size (0.4) / 2");
      desc.add<double>("dPtMaxFactor", 3)->setComment("dPt < 3 * resolution");
    }

    template <class T>
    void getTokens(const edm::Event& event) {
      if constexpr (std::is_same<T, Jet>::value) {
        // Skip token fetching conditionally for pat jets
        if (m_use_own_match) return;
      }
      event.getByToken(m_genJetsToken, m_genJets);
    }

    template <class T>
    const reco::GenJet* match(const T& jet, double resolution) {
      const reco::GenJet* matched_genJet = nullptr;

      if constexpr (std::is_same<T, Jet>::value) {
        if (m_use_own_match) {
          // The match is already there for pat jets, so we don't want to waste time looking.
          // This method is almost one-to-one inclusive for using slimmedGenJets.
          // Once in a few thousand events, there is a low-pt (~sub 15 GeV) jet that is not matched.
          // This is more of a philosophical question with the matching than a problem with the following code:
          // how come slimmedGenJets differs from the internal match?
          matched_genJet = jet.genJet();
          if (matched_genJet) {
            // The given conditions can be more harsh than those originally required for gen jets, so we check these.
            const double dPt = std::abs(matched_genJet->pt() - jet.pt());
            if (dPt > m_dPt_max_factor * resolution || deltaR(*matched_genJet, jet) >= m_dR_max) {
              matched_genJet = nullptr;
            }
          }
          // Return nullptr, if any of the conditions fails.
          return matched_genJet;
        }
      }
      // Only look at the gen jet collection if we need it.
      const reco::GenJetCollection& genJets = *m_genJets;

      // Try to find a gen jet matching
      // dR < m_dR_max
      // dPt < m_dPt_max_factor * resolution

      double min_dR = std::numeric_limits<double>::infinity();

      for (const auto& genJet : genJets) {
        const double dR = deltaR(genJet, jet);

        if (dR > min_dR)
          continue;

        if (dR < m_dR_max) {
          const double dPt = std::abs(genJet.pt() - jet.pt());
          if (dPt > m_dPt_max_factor * resolution)
            continue;

          min_dR = dR;
          matched_genJet = &genJet;
        }
      }

      return matched_genJet;
    }

  private:
    edm::InputTag m_genJetsTag;
    edm::EDGetTokenT<reco::GenJetCollection> m_genJetsToken;
    edm::Handle<reco::GenJetCollection> m_genJets;

    const bool m_use_own_match;
    const double m_dR_max;
    const double m_dPt_max_factor;
  };
};  // namespace pat

template <typename T>
class SmearedJetProducerT : public edm::stream::EDProducer<> {
  using JetCollection = std::vector<T>;

public:
  explicit SmearedJetProducerT(const edm::ParameterSet& cfg)
      : m_enabled(cfg.getParameter<bool>("enabled")),
        m_useDeterministicSeed(cfg.getParameter<bool>("useDeterministicSeed")),
        m_debug(cfg.getUntrackedParameter<bool>("debug", false)) {
    m_jets_token = consumes<JetCollection>(cfg.getParameter<edm::InputTag>("src"));

    if (m_enabled) {
      m_rho_token = consumes<double>(cfg.getParameter<edm::InputTag>("rho"));

      m_use_txt_files = cfg.exists("resolutionFile") && cfg.exists("scaleFactorFile");

      if (m_use_txt_files) {
        std::string resolutionFile = cfg.getParameter<edm::FileInPath>("resolutionFile").fullPath();
        std::string scaleFactorFile = cfg.getParameter<edm::FileInPath>("scaleFactorFile").fullPath();

        m_resolution_from_file.reset(new JME::JetResolution(resolutionFile));
        m_scale_factor_from_file.reset(new JME::JetResolutionScaleFactor(scaleFactorFile));
      } else {
        m_jets_algo = cfg.getParameter<std::string>("algo");
        m_jets_algo_pt = cfg.getParameter<std::string>("algopt");
      }

      std::uint32_t seed = cfg.getParameter<std::uint32_t>("seed");
      m_random_generator = std::mt19937(seed);

      bool skipGenMatching = cfg.getParameter<bool>("skipGenMatching");
      if (!skipGenMatching)
        m_genJetMatcher = std::make_shared<pat::GenJetMatcher>(cfg, consumesCollector());

      std::int32_t variation = cfg.getParameter<std::int32_t>("variation");
      m_uncertaintySource = cfg.getParameter<std::string>("uncertaintySource");
      m_nomVar = 1;
      if (variation == 0)
        m_systematic_variation = Variation::NOMINAL;
      else if (variation == 1)
        m_systematic_variation = Variation::UP;
      else if (variation == -1)
        m_systematic_variation = Variation::DOWN;
      else if (variation == 101) {
        m_systematic_variation = Variation::NOMINAL;
        m_nomVar = 1;
      } else if (variation == -101) {
        m_systematic_variation = Variation::NOMINAL;
        m_nomVar = -1;
      } else
        throw edm::Exception(edm::errors::ConfigFileReadError,
                             "Invalid value for 'variation' parameter. Only -1, 0, 1 or 101, -101 are supported.");
    }

    produces<JetCollection>();
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    // smearedJets
    edm::ParameterSetDescription desc;
    desc.setComment("A generic jet smearing module.");

    desc.add<edm::InputTag>("src")->setComment("Jet collection to smear.");
    desc.add<bool>("enabled", true)->setComment("If False, no smearing is performed.");
    desc.add<edm::InputTag>("rho")->setComment("Rho required for JER calculation.");
    desc.add<int>("variation", 0)->setComment("Systematic variation: 0: Nominal, -1: -sigma (down), 1: +sigma (up).");
    desc.add<std::string>("uncertaintySource", "")->setComment("If not specified, default to Total.");
    desc.add<bool>("useDeterministicSeed", true)->setComment("Pseudo-random seed based on event indices & such.");
    desc.add<unsigned int>("seed", 37428479)->setComment("Default seed, if useDeterministicSeed is set to false.");

    desc.add<bool>("skipGenMatching", false)->setComment("If True, always skip gen jet matching and smear jet with a random gaussian.");
    desc.addUntracked<bool>("debug", false);

    // The user must decide, which pair to use.
    auto source = (edm::ParameterDescription<std::string>("algo", true) and
                   edm::ParameterDescription<std::string>("algopt", true)) xor
                  (edm::ParameterDescription<edm::FileInPath>("resolutionFile", true) and
                   edm::ParameterDescription<edm::FileInPath>("scaleFactorFile", true));
    source->setComment("Read from GT: algo & algopt. From text files: resolutionFile & scaleFactorFile.");
    desc.addNode(std::move(source));

    pat::GenJetMatcher::fillDescriptions(desc);

    if constexpr (std::is_same<T, pat::Jet>::value) {
      descriptions.add("smearedPATJetProducer", desc);
    } else if constexpr (std::is_same<T, reco::PFJet>::value) {
      descriptions.add("smearedPFJetProducer", desc);
    } else if constexpr (std::is_same<T, reco::CaloJet>::value) {
      descriptions.add("smearedCaloJetProducer", desc);
    } else {
      // If we end up here, a new jet type has been added and the smearing code requires revision. 
      throw cms::Exception("BuildError") << "Unknown jet type for SmearedJetProducer!" << endl;
    }
  }

  void produce(edm::Event& event, const edm::EventSetup& setup) override {
    edm::Handle<JetCollection> jets_collection;
    event.getByToken(m_jets_token, jets_collection);

    // Disable the module when running on real data
    if (m_enabled && event.isRealData()) {
      m_enabled = false;
      m_genJetMatcher.reset();
    }

    edm::Handle<double> rho;
    if (m_enabled)
      event.getByToken(m_rho_token, rho);

    JME::JetResolution resolution;
    JME::JetResolutionScaleFactor resolution_sf;

    const JetCollection& jets = *jets_collection;

    if (m_enabled) {
      if (m_use_txt_files) {
        resolution = *m_resolution_from_file;
        resolution_sf = *m_scale_factor_from_file;
      } else {
        resolution = JME::JetResolution::get(setup, m_jets_algo_pt);
        resolution_sf = JME::JetResolutionScaleFactor::get(setup, m_jets_algo);
      }

      if (m_useDeterministicSeed) {
        unsigned int runNum_uint = static_cast<unsigned int>(event.id().run());
        unsigned int lumiNum_uint = static_cast<unsigned int>(event.id().luminosityBlock());
        unsigned int evNum_uint = static_cast<unsigned int>(event.id().event());
        unsigned int jet0eta = uint32_t(jets.empty() ? 0 : jets[0].eta() / 0.01);
        std::uint32_t seed = jet0eta + m_nomVar + (lumiNum_uint << 10) + (runNum_uint << 20) + evNum_uint;
        m_random_generator.seed(seed);
      }
    }

    if (m_genJetMatcher)
      m_genJetMatcher->getTokens<T>(event);

    auto smearedJets = std::make_unique<JetCollection>();

    for (const auto& jet : jets) {
      if ((!m_enabled) || (jet.pt() == 0)) {
        // Module disabled or invalid p4. Simply copy the input jet.
        smearedJets->push_back(jet);

        continue;
      }

      double jet_resolution = resolution.getResolution(
          {{JME::Binning::JetPt, jet.pt()}, {JME::Binning::JetEta, jet.eta()}, {JME::Binning::Rho, *rho}});
      double jer_sf = resolution_sf.getScaleFactor({{JME::Binning::JetPt, jet.pt()}, {JME::Binning::JetEta, jet.eta()}},
                                                   m_systematic_variation,
                                                   m_uncertaintySource);
      if (m_debug) {
        std::cout << "jet:  pt: " << jet.pt() << "  eta: " << jet.eta() << "  phi: " << jet.phi()
                  << "  e: " << jet.energy() << std::endl;
        std::cout << "resolution: " << jet_resolution << std::endl;
        std::cout << "resolution scale factor: " << jer_sf << std::endl;
      }

      const reco::GenJet* genJet = nullptr;
      if (m_genJetMatcher)
        genJet = m_genJetMatcher->match(jet, jet.pt() * jet_resolution);

      double smearFactor = 1.;

      if (genJet) {
        /*
                     * Case 1: we have a "good" gen jet matched to the reco jet
                     */

        if (m_debug) {
          std::cout << "gen jet:  pt: " << genJet->pt() << "  eta: " << genJet->eta() << "  phi: " << genJet->phi()
                    << "  e: " << genJet->energy() << std::endl;
        }

        double dPt = jet.pt() - genJet->pt();
        smearFactor = 1 + m_nomVar * (jer_sf - 1.) * dPt / jet.pt();
      } else if (jer_sf > 1) {
        /*
                     * Case 2: we don't have a gen jet. Smear jet pt using a random gaussian variation
                     */

        double sigma = jet_resolution * std::sqrt(jer_sf * jer_sf - 1);
        if (m_debug) {
          std::cout << "gaussian width: " << sigma << std::endl;
        }

        std::normal_distribution<> d(0, sigma);
        smearFactor = 1. + m_nomVar * d(m_random_generator);
      } else if (m_debug) {
        std::cout << "Impossible to smear this jet" << std::endl;
      }

      if (jet.energy() * smearFactor < MIN_JET_ENERGY) {
        // Negative or too small smearFactor. We would change direction of the jet
        // and this is not what we want.
        // Recompute the smearing factor in order to have jet.energy() == MIN_JET_ENERGY
        double newSmearFactor = MIN_JET_ENERGY / jet.energy();
        if (m_debug) {
          std::cout << "The smearing factor (" << smearFactor << ") is either negative or too small. Fixing it to "
                    << newSmearFactor << " to avoid change of direction." << std::endl;
        }
        smearFactor = newSmearFactor;
      }

      T smearedJet = jet;
      smearedJet.scaleEnergy(smearFactor);

      if (m_debug) {
        std::cout << "smeared jet (" << smearFactor << "):  pt: " << smearedJet.pt() << "  eta: " << smearedJet.eta()
                  << "  phi: " << smearedJet.phi() << "  e: " << smearedJet.energy() << std::endl;
      }

      smearedJets->push_back(smearedJet);
    }

    // Sort jets by pt
    std::sort(smearedJets->begin(), smearedJets->end(), jetPtComparator);

    event.put(std::move(smearedJets));
  }

private:
  static constexpr const double MIN_JET_ENERGY = 1e-2;

  edm::EDGetTokenT<JetCollection> m_jets_token;
  edm::EDGetTokenT<double> m_rho_token;
  bool m_enabled;
  std::string m_jets_algo_pt;
  std::string m_jets_algo;
  Variation m_systematic_variation;
  std::string m_uncertaintySource;
  bool m_useDeterministicSeed;
  bool m_debug;
  std::shared_ptr<pat::GenJetMatcher> m_genJetMatcher;

  bool m_use_txt_files;
  std::unique_ptr<JME::JetResolution> m_resolution_from_file;
  std::unique_ptr<JME::JetResolutionScaleFactor> m_scale_factor_from_file;

  std::mt19937 m_random_generator;

  GreaterByPt<T> jetPtComparator;

  int m_nomVar;
};
#endif
