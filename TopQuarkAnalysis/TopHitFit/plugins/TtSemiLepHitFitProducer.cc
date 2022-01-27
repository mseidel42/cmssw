#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "PhysicsTools/JetMCUtils/interface/combination.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiLepEvtPartons.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiEvtSolution.h"

#include "TopQuarkAnalysis/TopHitFit/interface/RunHitFit.h"
#include "TopQuarkAnalysis/TopHitFit/interface/LeptonTranslatorBase.h"
#include "TopQuarkAnalysis/TopHitFit/interface/JetTranslatorBase.h"
#include "TopQuarkAnalysis/TopHitFit/interface/METTranslatorBase.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

template <typename Lepton>
class TtSemiLepHitFitProducer : public edm::EDProducer {
  typedef edm::View<pat::Jet> JetCollection;
  typedef edm::View<pat::MET> METCollection;
  typedef edm::View<Lepton> LepCollection;
  typedef vector<pat::Particle> ParticleCollection;

public:
  explicit TtSemiLepHitFitProducer(const edm::ParameterSet &);
  ~TtSemiLepHitFitProducer() override;

private:
  // produce
  void produce(edm::Event &, const edm::EventSetup &) override;

  edm::EDGetTokenT<JetCollection> jetsToken_;
  edm::EDGetTokenT<METCollection> metsToken_;
  edm::EDGetTokenT<LepCollection> lepsToken_;

  struct FitResult {
    int Status;
    double Chi2;
    double Prob;
    double MT;
    double SigMT;
    pat::Particle HadB;
    pat::Particle HadP;
    pat::Particle HadQ;
    pat::Particle LepB;
    pat::Particle LepL;
    pat::Particle LepN;
    vector<int> JetCombi;
    int drop1;
    int drop2;
    bool operator<(const FitResult &rhs) { return Chi2 < rhs.Chi2; };
  };

  int runHitFit(std::list<FitResult> &FitResultList,
                edm::Handle<JetCollection> &jets,
                edm::Handle<METCollection> &mets,
                edm::Handle<LepCollection> &leps,
                size_t skipIdx = 4,
                size_t skipIdx2 = 5);

  typedef hitfit::RunHitFit<pat::Electron, pat::Muon, pat::Jet, pat::MET> PatHitFit;

  edm::FileInPath hfDefault_;
  edm::FileInPath hfElectronResolution_;
  edm::FileInPath hfMuonResolution_;
  edm::FileInPath hfUdscJetResolution_;
  edm::FileInPath hfBJetResolution_;
  edm::FileInPath hfMETResolution_;

  /// maximal number of jets (-1 possible to indicate 'all')
  const int maxNJets_;
  /// maximal number of combinations to be written to the event
  const int maxNComb_;

  /// input tag for b-tagging algorithm
  const string bTagAlgo_;
  /// min value of bTag for a b-jet
  const double minBTagValueBJet_;
  /// max value of bTag for a non-b-jet
  const double maxBTagValueLJet_;
  /// switch to tell whether to use b-tagging or not
  const bool useBTag_;
  /// switch to check the leading six jets instead of four
  const bool sixJets_;
  /// switch to check the leading five jets instead of four
  const bool fiveJets_;
  /// switch for selecting the two most likely b-jets
  const bool bestBs_;

  /// constraints
  const double mW_;
  const double mTop_;

  /// jet correction level
  const string jetCorrectionLevel_;

  bool firstRun_;

  const double probEpsilon_ = 1e-10;

  /*
    Jet permutation according to TQAF convention
    11 : leptonic b
    12 : hadronic b
    13 : hadronic W
    14 : hadronic W
  */
  const std::map<int, const int> type2idx_ = {{11, TtSemiLepEvtPartons::LepB},
                                              {12, TtSemiLepEvtPartons::HadB},
                                              {13, TtSemiLepEvtPartons::LightQ},
                                              {14, TtSemiLepEvtPartons::LightQBar}};

  hitfit::LeptonTranslatorBase<pat::Electron> electronTranslator_;
  hitfit::LeptonTranslatorBase<pat::Muon> muonTranslator_;
  hitfit::JetTranslatorBase<pat::Jet> jetTranslator_;
  hitfit::METTranslatorBase<pat::MET> metTranslator_;

  PatHitFit *HitFit;

public:
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
    // hitFitTtSemiLepEvent
    edm::ParameterSetDescription desc;

    desc.add<edm::InputTag>("jets", edm::InputTag("selectedPatJets"));
    desc.add<edm::InputTag>("mets", edm::InputTag("patMETs"));
    // ------------------------------------------------
    // maximum number of jets to be considered in the jet combinatorics
    // (has to be >= 4, can be set to -1 if you want to take all)
    // ------------------------------------------------
    desc.add<int>("maxNJets", 4);
    // ------------------------------------------------
    // maximum number of jet combinations finally written into the event, starting from the "best"
    // (has to be >= 1, can be set to -1 if you want to take all)
    // ------------------------------------------------
    desc.add<int>("maxNComb", 1);
    // ------------------------------------------------
    // option to use b-tagging
    // ------------------------------------------------
    desc.add<string>("bTagAlgo", "")->setComment("The user should absolutely provide this.");
    // Usually the following two are the same, but we allow them to be different
    desc.add<double>("minBDiscBJets", 0.0);
    desc.add<double>("maxBDiscLightJets", 1.0);
    // Important boolean flags
    desc.add<bool>("useBTagging", true);
    desc.add<bool>("useFiveJets", false);
    desc.add<bool>("useSixJets", false);
    desc.add<bool>("useBTagEmulation", false);
    // ------------------------------------------------
    // set mass values used in the constraints
    // set mass to 0 for no constraint
    // set mass to -1 to treat solve_nu_tmass and equal_side as false
    // ------------------------------------------------
    desc.add<double>("mW", 80.4);
    desc.add<double>("mTop", 0.0);
    // ------------------------------------------------
    // specify jet correction level as, Uncorrected, L1Offset, L2Relative, L3Absolute, L4Emf,
    // L5Hadron, L6UE, L7Parton. This is currently only used for debugging. The jet collection
    // should be "ready" when it arrives here.
    // ------------------------------------------------
    desc.add<string>("jetCorrectionLevel", "L3Absolute");
    // ------------------------------------------------
    // check that the jet collection behaves as expected, including e.g. the JEC level.
    // ------------------------------------------------
    desc.add<bool>("runJetTests", true);
    // ------------------------------------------------
    // resolutions and other configs
    // ------------------------------------------------
    desc.add<edm::FileInPath>("hitfitDefault",
                              edm::FileInPath("TopQuarkAnalysis/TopHitFit/data/setting/RunHitFitConfiguration.txt"));
    desc.add<edm::FileInPath>("hitfitElectronResolution",
                              edm::FileInPath("TopQuarkAnalysis/TopHitFit/data/resolution/tqafElectronResolution.txt"));
    desc.add<edm::FileInPath>("hitfitMuonResolution",
                              edm::FileInPath("TopQuarkAnalysis/TopHitFit/data/resolution/tqafMuonResolution.txt"));
    desc.add<edm::FileInPath>("hitfitUdscJetResolution",
                              edm::FileInPath("TopQuarkAnalysis/TopHitFit/data/resolution/tqafUdscJetResolution.txt"));
    desc.add<edm::FileInPath>("hitfitBJetResolution",
                              edm::FileInPath("TopQuarkAnalysis/TopHitFit/data/resolution/tqafBJetResolution.txt"));
    desc.add<edm::FileInPath>("hitfitMETResolution",
                              edm::FileInPath("TopQuarkAnalysis/TopHitFit/data/resolution/tqafKtResolution.txt"));

    // Only electrons and muons are catered
    if constexpr (std::is_same<Lepton, pat::Electron>::value) {
      desc.add<edm::InputTag>("leps", edm::InputTag("selectedPatElectrons"));
      descriptions.add("hitFitTtSemiLepEventElectrons", desc);
    }
    if constexpr (std::is_same<Lepton, pat::Muon>::value) {
      desc.add<edm::InputTag>("leps", edm::InputTag("selectedPatMuons"));
      descriptions.add("hitFitTtSemiLepEventMuons", desc);
    }
  }
};

template <typename Lepton>
TtSemiLepHitFitProducer<Lepton>::TtSemiLepHitFitProducer(const edm::ParameterSet &cfg)
    : jetsToken_(consumes<JetCollection>(cfg.getParameter<edm::InputTag>("jets"))),
      metsToken_(consumes<METCollection>(cfg.getParameter<edm::InputTag>("mets"))),
      lepsToken_(consumes<LepCollection>(cfg.getParameter<edm::InputTag>("leps"))),

      // The following five initializers read the config parameters for the: ASCII text files which contains the physics object resolutions.
      hfDefault_(cfg.getParameter<edm::FileInPath>("hitfitDefault")),
      hfElectronResolution_(cfg.getParameter<edm::FileInPath>("hitfitElectronResolution")),
      hfMuonResolution_(cfg.getParameter<edm::FileInPath>("hitfitMuonResolution")),
      hfUdscJetResolution_(cfg.getParameter<edm::FileInPath>("hitfitUdscJetResolution")),
      hfBJetResolution_(cfg.getParameter<edm::FileInPath>("hitfitBJetResolution")),
      hfMETResolution_(cfg.getParameter<edm::FileInPath>("hitfitMETResolution")),

      // Constants
      maxNJets_(cfg.getParameter<int>("maxNJets")),
      maxNComb_(cfg.getParameter<int>("maxNComb")),
      bTagAlgo_(cfg.getParameter<string>("bTagAlgo")),
      minBTagValueBJet_(cfg.getParameter<double>("minBDiscBJets")),
      maxBTagValueLJet_(cfg.getParameter<double>("maxBDiscLightJets")),
      useBTag_(cfg.getParameter<bool>("useBTagging") && bTagAlgo_ != ""),
      sixJets_(cfg.getParameter<bool>("useSixJets")),
      fiveJets_(cfg.getParameter<bool>("useFiveJets") || sixJets_),
      bestBs_(cfg.getParameter<bool>("useBTagEmulation")),
      mW_(cfg.getParameter<double>("mW")),
      mTop_(cfg.getParameter<double>("mTop")),
      jetCorrectionLevel_(cfg.getParameter<string>("jetCorrectionLevel")),
      firstRun_(cfg.getParameter<bool>("runJetTests")),

      // The following four initializers instantiate the translator between PAT objects and HitFit objects using the ASCII text files which contains the resolutions.
      electronTranslator_(hfElectronResolution_.fullPath()),
      muonTranslator_(hfMuonResolution_.fullPath()),
      jetTranslator_(hfUdscJetResolution_.fullPath(), hfBJetResolution_.fullPath()),
      metTranslator_(hfMETResolution_.fullPath()) {
  // Create an instance of RunHitFit and initialize it.
  HitFit = new PatHitFit(electronTranslator_,
                         muonTranslator_,
                         jetTranslator_,
                         metTranslator_,
                         hfDefault_.fullPath(),
                         mW_,
                         mW_,
                         mTop_,
                         useBTag_);

  edm::LogVerbatim("TopHitFit") << "\n"
                                << "+++++++++++ TtSemiLepHitFitProducer ++++++++++++ \n"
                                << " Due to the eta ranges for which resolutions     \n"
                                << " are provided in                                 \n"
                                << " TopQuarkAnalysis/TopHitFit/data/resolution/     \n"
                                << " so far, the following cuts are currently        \n"
                                << " implemented in the TtSemiLepHitFitProducer:     \n"
                                << "++++++++++++++++++++++++++++++++++++++++++++++++ \n";

  // Informatics!
  cout << endl << endl;
  if (useBTag_) {
    if (bestBs_) {
      cout << "||| Running HitFit in b-tag emulation mode! |||" << endl;
    } else {
      cout << "||| Running HitFit in b-tagging mode with limits " << maxBTagValueLJet_ << "/" << minBTagValueBJet_
           << " |||" << endl;
      if (minBTagValueBJet_ < maxBTagValueLJet_) {
        throw cms::Exception("TtSemiLepHitFitProducer") << "ERROR!!! Overlapping tags!!!";
      }
    }
  } else {
    cout << "||| Running HitFit in brute force mode, without b-tagging! |||" << endl;
  }
  if (fiveJets_) {
    if (sixJets_) {
      cout << "||| Permuting through the six leading jets instead of four! |||" << endl;
    } else {
      cout << "||| Permuting through the five leading jets instead of four! |||" << endl;
    }
  } else {
    cout << "||| Permuting through the four leading jets! |||" << endl;
  }
  cout << "||| Probability epsilon for throwing out sub-leading permutations within the chosen " << maxNComb_
       << " permutations: " << probEpsilon_ << " |||" << endl;
  cout << endl << endl;

  produces<ParticleCollection>("PartonsHadP");
  produces<ParticleCollection>("PartonsHadQ");
  produces<ParticleCollection>("PartonsHadB");
  produces<ParticleCollection>("PartonsLepB");
  produces<ParticleCollection>("Leptons");
  produces<ParticleCollection>("Neutrinos");

  produces<vector<vector<int> > >();
  produces<vector<double> >("Chi2");
  produces<vector<double> >("Prob");
  produces<vector<double> >("MT");
  produces<vector<double> >("SigMT");
  produces<vector<int> >("Status");
  produces<int>("NumberOfConsideredJets");
  produces<int>("Drop1");
  produces<int>("Drop2");
}

template <typename Lepton>
TtSemiLepHitFitProducer<Lepton>::~TtSemiLepHitFitProducer() {
  delete HitFit;
}

template <typename Lepton>
int TtSemiLepHitFitProducer<Lepton>::runHitFit(std::list<FitResult> &FitResultList,
                                               edm::Handle<JetCollection> &jets,
                                               edm::Handle<METCollection> &mets,
                                               edm::Handle<LepCollection> &leps,
                                               size_t skipIdx,
                                               size_t skipIdx2) {
  if (useBTag_) {
    // In this running mode, the two leading jets are reserved for b-jets, and should not be skipped.
    if (skipIdx <= 1) {
      throw cms::Exception("TtSemiLepHitFitProducer") << "The B-jets should not be skipped!";
    }
  }
  if (skipIdx >= skipIdx2) {
    throw cms::Exception("TtSemiLepHitFitProducer") << "Wrong order for skipping indices!";
  }
  // Clear the internal state
  HitFit->clear();

  // Add lepton into HitFit
  bool foundLepton = false;
  for (size_t iLep = 0; iLep < leps->size(); ++iLep) {
    // Eta sanity check: the leptons should already be sanitized.
    if (fabs(leps->at(iLep).eta()) <= 2.5) {
      HitFit->AddLepton((*leps)[iLep]);
      foundLepton = true;
    } else {
      // Fatal issue: this should not occur, so we make some noise.
      throw cms::Exception("TtSemiLepHitFitProducer")
          << "Issues with lepton eta. " << fabs(leps->at(iLep).eta()) << " vs. max. 2.5";
    }
    break;
  }
  if (!foundLepton)
    return 0;

  // Add jets into HitFit
  int nJetsFound = 0, nBJetsFound = 0, nLJetsFound = 0;
  for (size_t iJet = 0; iJet < jets->size(); ++iJet) {
    // Allow skipping a certain jets
    if (iJet == skipIdx || iJet == skipIdx2)
      continue;
    const auto &jet = jets->at(iJet);
    const double jet_abseta = std::abs(jet.eta());
    // If Calo jet usage is resurrected: jet.isCaloJet() ? ((reco::CaloJet*) jet.originalObject())->detectorP4().eta()
    // Eta sanity check: the jets should already be sanitized
    if (jet_abseta <= 5.2) {
      if (useBTag_) {
        const double bTag = jet.bDiscriminator(bTagAlgo_);
        // In the bestBs_ mode, the two leading jets are the best b-jet candidates.
        const bool isB = bestBs_ ? (nJetsFound < 2) : bTag > minBTagValueBJet_;
        const bool isL = bestBs_ ? (!isB) : bTag <= maxBTagValueLJet_;
        HitFit->AddJet(jet, isB, isL);
        if (isB)
          ++nBJetsFound;
        if (isL)
          ++nLJetsFound;
      } else {
        // With no b-tagging, we are less elaborate and use more brute force to loop over permutations.
        HitFit->AddJet(jet);
      }
      if (++nJetsFound == maxNJets_)
        break;
    } else {
      // Fatal issue: this should not occur, so we make some noise.
      throw cms::Exception("TtSemiLepHitFitProducer") << "Issues with jet eta. " << jet_abseta << " v. max. 5.2";
    }
  }

  if (nJetsFound < 4 || (useBTag_ && (nBJetsFound < 2 || nLJetsFound < 2)))
    return nJetsFound;

  // Add missing transverse energy into HitFit
  HitFit->SetMet((*mets)[0]);

  //
  // R U N   H I T F I T
  //
  // Run the kinematic fit and get how many permutations are possible
  // in the fit
  HitFit->FitAllPermutation();

  // Get the fit results for all permutations
  vector<hitfit::Fit_Result> hitFitResult = HitFit->GetFitAllPermutation();

  //
  // BEGIN PART WHICH EXTRACTS INFORMATION FROM HITFIT
  //

  // Loop over all permutations and extract the information, save into a vector that is later sorted according to chi2
  for (const auto &fitProd : hitFitResult) {
    const auto &fit = fitProd.ev();
    vector<int> hitCombi(4, -1);
    if (fit.njets() > 4) {
      // Fatal issue: this should not occur, so we make some noise.
      throw cms::Exception("TtSemiLepHitFitProducer") << "Too many jets in the fit!";
    }

    // Get the number of jets and loop over the jets
    for (size_t jet = 0, nJets = fit.njets(); jet < nJets; ++jet) {
      int jet_type = fit.jet(jet).type();

      if (jet_type > 10 && jet_type < 15)
        hitCombi[type2idx_.at(jet_type)] = jet;
    }

    // Store the kinematic quantities in the corresponding containers.
    hitfit::Lepjets_Event_Jet lepB = fit.jet(hitCombi[type2idx_.at(11)]);
    hitfit::Lepjets_Event_Jet hadB = fit.jet(hitCombi[type2idx_.at(12)]);
    hitfit::Lepjets_Event_Jet hadP = fit.jet(hitCombi[type2idx_.at(13)]);
    hitfit::Lepjets_Event_Jet hadQ = fit.jet(hitCombi[type2idx_.at(14)]);
    hitfit::Lepjets_Event_Lep lepL = fit.lep(0);
    const auto &metL = fit.met();
    // Later on, hitCombi will be used for the full jet collection, so correcting for this at skipIdx.
    if (skipIdx < 4) {
      int skipInt = skipIdx;
      for (auto &hc : hitCombi) {
        if (hc >= skipInt)
          hc += 1;
      }
    }
    if (skipIdx2 < 5) {
      int skipInt = skipIdx2;
      for (auto &hc : hitCombi) {
        if (hc >= skipInt)
          hc += 1;
      }
    }

    if (fitProd.chisq() > 0) {
      FitResultList.emplace_back();
      auto &result = FitResultList.back();

      result.Status = 0;
      result.Chi2 = fitProd.chisq();
      result.Prob = exp(-1.0 * (fitProd.chisq()) / 2.0);
      result.MT = fitProd.mt();
      result.SigMT = fitProd.sigmt();
      result.HadB = pat::Particle(reco::LeafCandidate(
          0, math::XYZTLorentzVector(hadB.p().x(), hadB.p().y(), hadB.p().z(), hadB.p().t()), math::XYZPoint()));
      result.HadP = pat::Particle(reco::LeafCandidate(
          0, math::XYZTLorentzVector(hadP.p().x(), hadP.p().y(), hadP.p().z(), hadP.p().t()), math::XYZPoint()));
      result.HadQ = pat::Particle(reco::LeafCandidate(
          0, math::XYZTLorentzVector(hadQ.p().x(), hadQ.p().y(), hadQ.p().z(), hadQ.p().t()), math::XYZPoint()));
      result.LepB = pat::Particle(reco::LeafCandidate(
          0, math::XYZTLorentzVector(lepB.p().x(), lepB.p().y(), lepB.p().z(), lepB.p().t()), math::XYZPoint()));
      result.LepL = pat::Particle(reco::LeafCandidate(
          0, math::XYZTLorentzVector(lepL.p().x(), lepL.p().y(), lepL.p().z(), lepL.p().t()), math::XYZPoint()));
      result.LepN = pat::Particle(
          reco::LeafCandidate(0, math::XYZTLorentzVector(metL.x(), metL.y(), metL.z(), metL.t()), math::XYZPoint()));
      result.JetCombi = hitCombi;
      result.drop1 = skipIdx;
      result.drop2 = skipIdx2;
    }
  }
  return nJetsFound;
}

template <typename Lepton>
void TtSemiLepHitFitProducer<Lepton>::produce(edm::Event &evt, const edm::EventSetup &setup) {
  auto jets = evt.getHandle(jetsToken_);
  auto mets = evt.getHandle(metsToken_);
  auto leps = evt.getHandle(lepsToken_);

  // We run a bunch of tests for the delivered jets at the first run.
  if (firstRun_) {
    if (!jets.isValid()) {
      throw cms::Exception("TtSemiLepHitFitProducer") << "Faulty jet collection provided!";
    }
    if (!mets.isValid()) {
      throw cms::Exception("TtSemiLepHitFitProducer") << "Faulty MET collection provided!";
    }
    if (!leps.isValid()) {
      throw cms::Exception("TtSemiLepHitFitProducer") << "Faulty lepton collection provided!";
    }
    if (jets->size() > 0) {
      cout << endl << endl << "##################################" << endl;
      cout << "Testing the jets within HitFit... If any errors occur, start by checking at ";
      cout << "TopQuarkAnalysis/TopHitFit/plugins/TtSemiLepHitFitProducer.h!" << endl;
      for (size_t iJet = 0; iJet < jets->size(); ++iJet) {
        const auto &jet = jets->at(iJet);
        // Default: jet.isPFJet(), otherwise the user must edit this code to show she/he knows what they are doing.
        // Automatic swithcing is quite dangerous, as the user might not understand what they are doing.
        if (!jet.isPFJet()) {
          cout << "Caution: it is currently expected that PFJets are used in the kinematic fit!" << endl;
          if (jet.isCaloJet()) {
            // One needs to change the eta value from jet.eta() to static_cast<const reco::CaloJet*>(jet.originalObject())->detectorP4().eta().
            // This is done within this file, and in PatJetHitFitTranslator. The usage of CaloJets in a precision measurement is unexpected
            // currently, and the feature is disabled as it poses more a threat for bugs than any positive outcomes.
            throw cms::Exception("TtSemiLepHitFitProducer")
                << "It seems that Calo jets are used. Please visit"
                << "TopQuarkAnalysis/TopHitFit/src/PatJetHitFitTranslator.cc"
                << " and make the necessary edits, if you really want to proceed. Exiting!";
          } else {
            cout << "An unknown jet type was encountered! Is this PUPPI? Check if this affects HitFit." << endl;
            // One should probably check PatJetHitFitTranslator also in this case.
          }
        }
        // Check that the jet correction lavel is the same as given.
        if (jetCorrectionLevel_ != "") {
          auto sets = jet.availableJECSets();
          if (sets.size() < 1) {
            for (auto &set : sets)
              cout << set << " ";
            cout << endl;
            throw cms::Exception("TtSemiLepHitFitProducer")
                << "Maximally one JEC set expected, found " << sets.size() << "; please check!";
          } else if (sets.size() > 1) {
            cout << "Exotic: more than one JEC set found " << sets.size() << endl;
          }
          auto lvls = jet.availableJECLevels();
          if (lvls.size() < 2) {
            for (auto &lvl : lvls)
              cout << lvl << " ";
            cout << endl;
            throw cms::Exception("TtSemiLepHitFitProducer")
                << "Expected at least a couple levels of JEC, please check!";
          }
          if (jet.currentJECLevel() != jetCorrectionLevel_) {
            throw cms::Exception("TtSemiLepHitFitProducer")
                << "The JEC level in the jet collection " << jet.currentJECLevel()
                << " disagrees with the requested one, " << jetCorrectionLevel_ << ", please check! Exiting!";
          }
        }
      }
    }
    cout << "Jet tests performed on " << jets->size() << " jets OK!" << endl;
    cout << "##################################" << endl << endl;
    firstRun_ = false;
  }

  // -----------------------------------------------------
  // skip events with no appropriate lepton candidate in
  // or empty MET or less jets than partons
  // -----------------------------------------------------
  std::list<FitResult> FitResultList;
  int jetsConsidered = -1;
  // With b-tagging, each four jets produce 4 combinations and without, 24 combinations.
  // Base mode: 4/24 combinations.
  // Five jets: 12/120 combinations.
  // Six jets: 24/360 combinations.
  if (!jets->empty() && !mets->empty() && !leps->empty()) {
    // Jets 0, 1, 2, 3: 4/24 (b-tagging/not)
    jetsConsidered = runHitFit(FitResultList, jets, mets, leps);
    if (fiveJets_ && jets->size() > 4) {
      // We know now that there are at least 5 good jets
      jetsConsidered = 5;
      // Jets 0, 1, 2, 4: 4/24 (b-tagging/not)
      runHitFit(FitResultList, jets, mets, leps, 3);
      // Jets 0, 1, 3, 4: 4/24 (b-tagging/not)
      runHitFit(FitResultList, jets, mets, leps, 2);
      if (!useBTag_) {
        // Jets 0, 2, 3, 4: 4/24 (b-tagging/not)
        runHitFit(FitResultList, jets, mets, leps, 1);
        // Jets 1, 2, 3, 4: 4/24 (b-tagging/not)
        runHitFit(FitResultList, jets, mets, leps, 0);
      }
      if (sixJets_ && jets->size() > 5) {
        // We know now that there are at least 6 good jets
        jetsConsidered = 6;
        // Jets 0, 1, 2, 5: 4/24 (b-tagging/not)
        runHitFit(FitResultList, jets, mets, leps, 3, 4);
        // Jets 0, 1, 3, 4: 4/24 (b-tagging/not)
        runHitFit(FitResultList, jets, mets, leps, 2, 4);
        // Jets 0, 1, 4, 5: 4/24 (b-tagging/not)
        runHitFit(FitResultList, jets, mets, leps, 2, 3);
        if (!useBTag_) {
          // Jets 0, 2, 3, 5: 4/24 (b-tagging/not)
          runHitFit(FitResultList, jets, mets, leps, 1, 4);
          // Jets 0, 2, 4, 5: 4/24 (b-tagging/not)
          runHitFit(FitResultList, jets, mets, leps, 1, 3);
          // Jets 0, 3, 4, 5: 4/24 (b-tagging/not)
          runHitFit(FitResultList, jets, mets, leps, 1, 2);
          // Jets 1, 2, 3, 5: 4/24 (b-tagging/not)
          runHitFit(FitResultList, jets, mets, leps, 0, 4);
          // Jets 1, 2, 4, 5: 4/24 (b-tagging/not)
          runHitFit(FitResultList, jets, mets, leps, 0, 3);
          // Jets 1, 3, 4, 5: 4/24 (b-tagging/not)
          runHitFit(FitResultList, jets, mets, leps, 0, 2);
          // Jets 2, 3, 4, 5: 4/24 (b-tagging/not)
          runHitFit(FitResultList, jets, mets, leps, 0, 1);
        }
      }
    }
  }

  // -----------------------------------------------------
  // feed out result
  // starting with the JetComb having the smallest chi2
  // -----------------------------------------------------
  auto pPartonsHadP = std::make_unique<ParticleCollection>();
  auto pPartonsHadQ = std::make_unique<ParticleCollection>();
  auto pPartonsHadB = std::make_unique<ParticleCollection>();
  auto pPartonsLepB = std::make_unique<ParticleCollection>();
  auto pLeptons = std::make_unique<ParticleCollection>();
  auto pNeutrinos = std::make_unique<ParticleCollection>();

  auto pCombi = std::make_unique<vector<vector<int> > >();
  auto pChi2 = std::make_unique<vector<double> >();
  auto pProb = std::make_unique<vector<double> >();
  auto pMT = std::make_unique<vector<double> >();
  auto pSigMT = std::make_unique<vector<double> >();
  auto pStatus = std::make_unique<vector<int> >();

  int drop1 = -1;
  int drop2 = -1;
  if (FitResultList.size() == 0) {  // in case no fit results were stored in the list
    // the kinFit getters return empty objects here
    pPartonsHadP->emplace_back(pat::Particle());
    pPartonsHadQ->emplace_back(pat::Particle());
    pPartonsHadB->emplace_back(pat::Particle());
    pPartonsLepB->emplace_back(pat::Particle());
    pLeptons->emplace_back(pat::Particle());
    pNeutrinos->emplace_back(pat::Particle());
    // indices referring to the jet combination
    pCombi->emplace_back(4, -1);
    // chi2
    pChi2->push_back(-1.);
    // chi2 probability
    pProb->push_back(-1.);
    // fitted top mass
    pMT->push_back(-1.);
    pSigMT->push_back(-1.);
    // status of the fitter
    pStatus->push_back(-1);
  } else {
    // sort results w.r.t. chi2 values in increasing order
    FitResultList.sort();

    int iComb = 0;
    for (auto result = FitResultList.begin(); result != FitResultList.end(); ++result) {
      if (++iComb == 1) {
        drop1 = result->drop1;
        drop2 = result->drop2;
      } else if (result->Prob < probEpsilon_)
        break;
      if (maxNComb_ >= 0 and iComb > maxNComb_)
        break;
      // physics objects
      pPartonsHadP->emplace_back(result->HadP);
      pPartonsHadQ->emplace_back(result->HadQ);
      pPartonsHadB->emplace_back(result->HadB);
      pPartonsLepB->emplace_back(result->LepB);
      pLeptons->emplace_back(result->LepL);
      pNeutrinos->emplace_back(result->LepN);
      // indices referring to the jet combination
      pCombi->push_back(result->JetCombi);
      // chi2
      pChi2->push_back(result->Chi2);
      // chi2 probability
      pProb->push_back(result->Prob);
      // fitted top mass
      pMT->push_back(result->MT);
      pSigMT->push_back(result->SigMT);
      // status of the fitter
      pStatus->push_back(result->Status);
    }
  }
  evt.put(std::move(pCombi));
  evt.put(std::move(pPartonsHadP), "PartonsHadP");
  evt.put(std::move(pPartonsHadQ), "PartonsHadQ");
  evt.put(std::move(pPartonsHadB), "PartonsHadB");
  evt.put(std::move(pPartonsLepB), "PartonsLepB");
  evt.put(std::move(pLeptons), "Leptons");
  evt.put(std::move(pNeutrinos), "Neutrinos");
  evt.put(std::move(pChi2), "Chi2");
  evt.put(std::move(pProb), "Prob");
  evt.put(std::move(pMT), "MT");
  evt.put(std::move(pSigMT), "SigMT");
  evt.put(std::move(pStatus), "Status");
  evt.put(std::make_unique<int>(jetsConsidered), "NumberOfConsideredJets");
  evt.put(std::make_unique<int>(drop1), "Drop1");
  evt.put(std::make_unique<int>(drop2), "Drop2");
}

// Specifications!
typedef TtSemiLepHitFitProducer<pat::Muon> TtSemiLepHitFitProducerMuon;
typedef TtSemiLepHitFitProducer<pat::Electron> TtSemiLepHitFitProducerElectron;

DEFINE_FWK_MODULE(TtSemiLepHitFitProducerMuon);
DEFINE_FWK_MODULE(TtSemiLepHitFitProducerElectron);
