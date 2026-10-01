// Reads the user-facing logical truth::Graph (the GEN+SIM merged truth produced
// by TruthLogicalGraphProducer, the heart of the new "Truth Graph" feature) and
// rebuilds a HepMC3 GenEvent from it so that the same record can be fed to Rivet
// (RivetAnalyzer / ParticleLevelProducer) exactly as the standard
// GenParticles2HepMCConverter lets it be fed from the GEN/SIM reco::GenParticles.
//
// This is the truth-graph twin of the standard GenParticles2HepMCConverter.cc.
// Where that one reads reco::GenParticles + GenEventInfoProduct and re-erects the
// HepMC3 record, this one reads a truth::Graph + GenEventInfoProduct and does the
// same. The graph already merged the GEN and SIM sides and kept the standalone
// four-momentum and position for every node, so the conversion is a faithful walk
// of the bipartite Particle <-> Vertex graph.
//
// What is emitted by default mirrors the scope of the GEN record (the same scope
// the standard converter emits): the particles that carry a GEN back-reference,
// and the (non-artificial) vertices with a GEN back-reference. SIM-only
// particles/vertices are skipped because their four-momentum is a Geant4 quantity,
// not a generator one, which is not what Rivet analyses consume. The optional
// includeSimOnlyParticles flag lets the user add them anyway (e.g. for a
// truth-particle-level study that wants the Geant4 secondaries).
//
// The truth graph stores positions in (cm, ns); the HepMC3 default units are
// (mm, mm-of-c*t). Lengths are scaled 10x and the time is converted ns -> mm-of-c
// (299.792458 mm/ns), so the generator's vertex time is preserved - an improvement
// over the standard converter, which had to set t=0 because the reco::GenParticle
// does not carry one.

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/TruthInfo/interface/Graph.h"
#include "SimDataFormats/TruthInfo/interface/ParticleData.h"
#include "SimDataFormats/TruthInfo/interface/VertexData.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMC3Product.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Units.h"
#include "HepMC3/Print.h"
#include "HepMC3/WriterAscii.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <unordered_map>
#include <vector>

using namespace std;

class TruthGraph2HepMCConverter : public edm::stream::EDProducer<edm::stream::WatchRuns> {
public:
  explicit TruthGraph2HepMCConverter(const edm::ParameterSet& pset);
  ~TruthGraph2HepMCConverter() override;

  void beginRun(edm::Run const& iRun, edm::EventSetup const&) override;
  void produce(edm::Event& event, const edm::EventSetup& eventSetup) override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // Decision: should this truth-graph particle be emitted into the HepMC3 record?
  bool keepParticle(const truth::ParticleData& p) const;
  // Decision: should this truth-graph vertex be emitted into the HepMC3 record?
  bool keepVertex(const truth::VertexData& v) const;

  // Convert the graph's (cm, ns) position into a HepMC3 FourVector in (mm, mm-of-ct),
  // matching the standard convention of the rest of the code base.
  inline HepMC3::FourVector fourVector(const math::XYZTLorentzVectorD& p4_or_pos, bool isPosition) const {
    if (!isPosition) {
      // A particle four-momentum: the graph already stores GeV, which is the HepMC3 default.
      return HepMC3::FourVector(p4_or_pos.px(), p4_or_pos.py(), p4_or_pos.pz(), p4_or_pos.e());
    }
    // A vertex position: graph stores (cm, ns); HepMC3 default is (mm, mm-of-c*t).
    constexpr double cmToMm = 10.0;
    constexpr double nsToMmOfCT = 299.792458;  // 1 ns * c, in mm
    return HepMC3::FourVector(p4_or_pos.x() * cmToMm,
                              p4_or_pos.y() * cmToMm,
                              p4_or_pos.z() * cmToMm,
                              p4_or_pos.t() * nsToMmOfCT);
  }

  edm::EDGetTokenT<truth::Graph> graphToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
  edm::EDGetTokenT<GenRunInfoProduct> genRunInfoToken_;  // consumed InRun
  edm::ESGetToken<HepPDT::ParticleDataTable, PDTRecord> pTable_;

  const bool includeSimOnlyParticles_;  // emit SIM-only particles too (off by default)
  const bool includeArtificialVertices_;  // emit InitialState/UnderlyingEvent/... vertices too
  const double cmEnergy_;                  // fallback beam energy (when the graph has no beam protons)
  const bool writeHepMC_;                  // write a HepMC3 text file
  const std::string outputFile_;           // path of the HepMC3 text file
  const unsigned verbosity_;              // 0 silent, 1: a one-line summary per event

  HepMC3::GenCrossSectionPtr xsec_;
  std::shared_ptr<HepMC3::Writer> writer_;
};

TruthGraph2HepMCConverter::TruthGraph2HepMCConverter(const edm::ParameterSet& pset)
    : includeSimOnlyParticles_(pset.getUntrackedParameter<bool>("includeSimOnlyParticles", false)),
      includeArtificialVertices_(pset.getUntrackedParameter<bool>("includeArtificialVertices", false)),
      cmEnergy_(pset.getUntrackedParameter<double>("cmEnergy", 14000.)),
      writeHepMC_(pset.getUntrackedParameter<bool>("writeHepMC", false)),
      outputFile_(pset.getUntrackedParameter<std::string>("outputFile", "events.hepmc")),
      verbosity_(pset.getUntrackedParameter<unsigned>("verbosity", 0)) {
  graphToken_ = consumes<truth::Graph>(pset.getParameter<edm::InputTag>("truthGraph"));
  genEventInfoToken_ = consumes<GenEventInfoProduct>(pset.getParameter<edm::InputTag>("genEventInfo"));
  genRunInfoToken_ = consumes<GenRunInfoProduct, edm::InRun>(pset.getParameter<edm::InputTag>("genEventInfo"));
  pTable_ = esConsumes<HepPDT::ParticleDataTable, PDTRecord>();

  produces<edm::HepMC3Product>("unsmeared");

  if (writeHepMC_) {
    writer_ = std::make_shared<HepMC3::WriterAscii>(outputFile_);
  }
}

TruthGraph2HepMCConverter::~TruthGraph2HepMCConverter() {
  if (writeHepMC_) {
    writer_->close();
    writer_.reset();
  }
}

void TruthGraph2HepMCConverter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("truthGraph", edm::InputTag("truthLogicalGraphProducer"))
      ->setComment("the truth::Graph product (typically the output of TruthLogicalGraphProducer)");
  desc.add<edm::InputTag>("genEventInfo", edm::InputTag("generator"))
      ->setComment(
          "generator module label providing the GenEventInfoProduct (signalProcessID, qScale, "
          "alphaQCD, alphaQED, weights, pdf info) AND the GenRunInfoProduct (cross section).");

  desc.addUntracked<bool>("includeSimOnlyParticles", false)
      ->setComment(
          "Also emit SIM-only logical particles (Geant4 secondaries) in the HepMC3 record. Off by "
          "default, since their four-momentum is a Geant4 quantity and not a generator one, which "
          "is not what Rivet analyses consume. The GEN-side of a merged GEN+SIM particle is always "
          "emitted, so the standard GEN record is faithfully reproduced.");
  desc.addUntracked<bool>("includeArtificialVertices", false)
      ->setComment(
          "Also emit the artificial InitialState / UnderlyingEvent / Interaction / "
          "BeamSideInput vertices the selection post-processor may add. Off by default: those "
          "vertices are a bookkeeping device, not generator vertices, and Rivet expects a real "
          "HepMC3 record. Useful for debugging when the graph was built from a heavily pruned "
          "selection preset.");
  desc.addUntracked<double>("cmEnergy", 14000.)
      ->setComment(
          "Centre-of-mass energy (GeV) used as a fallback to synthesize the two beam protons "
          "when the input graph does not carry them (e.g. a particle-gun sample, or a graph "
          "produced with a selection preset that pruned the beam protons).");
  desc.addUntracked<bool>("writeHepMC", false)
      ->setComment("if true, also write the produced HepMC3 events to outputFile (for debugging).");
  desc.addUntracked<std::string>("outputFile", "events.hepmc")
      ->setComment("path of the HepMC3 text file when writeHepMC is true.");
  desc.addUntracked<unsigned>("verbosity", 0)
      ->setComment("0 = silent, >0 = print a one-line summary per event.");
  descriptions.addWithDefaultLabel(desc);
}

void TruthGraph2HepMCConverter::beginRun(edm::Run const& iRun, edm::EventSetup const&) {
  edm::Handle<GenRunInfoProduct> genRunInfoHandle;
  iRun.getByToken(genRunInfoToken_, genRunInfoHandle);

  xsec_ = make_shared<HepMC3::GenCrossSection>();
  if (genRunInfoHandle.isValid()) {
    xsec_->set_cross_section(genRunInfoHandle->internalXSec().value(), genRunInfoHandle->internalXSec().error());
  } else {
    // dummy cross section
    xsec_->set_cross_section(1., 0.);
  }
}

bool TruthGraph2HepMCConverter::keepParticle(const truth::ParticleData& p) const {
  // Synthetic particles (a connector the selection adds to attach a root to an
  // artificial Interaction vertex, or a SignalStandIn standing in for a resonance
  // the generator never wrote) are not generator particles: their momentum is an
  // accounting sum, not a generator quantity. Never emit them.
  if (p.isSynthetic()) {
    return false;
  }

  // By default only emit the GEN-side of the record, mirroring the scope of the
  // standard GenParticles2HepMCConverter (which fed Rivet from the GEN reco
  // particles). The GEN-side of a merged GEN+SIM particle is what a Rivet analysis
  // needs; the SIM-only Geant4 secondaries are not a generator quantity.
  if (!p.hasGen()) {
    return includeSimOnlyParticles_;
  }

  return true;
}

bool TruthGraph2HepMCConverter::keepVertex(const truth::VertexData& v) const {
  // Artificial vertices (InitialState / UnderlyingEvent / Interaction / BeamSideInput)
  // are the post-processing bookkeeping, not generator vertices. Skip them.
  if (v.isArtificial()) {
    return includeArtificialVertices_;
  }

  if (!v.hasGen()) {
    return includeSimOnlyParticles_;
  }

  return true;
}

void TruthGraph2HepMCConverter::produce(edm::Event& event, const edm::EventSetup& eventSetup) {
  auto const& graph = event.get(graphToken_);

  edm::Handle<GenEventInfoProduct> genEventInfoHandle;
  event.getByToken(genEventInfoToken_, genEventInfoHandle);

  auto const& pTableData = eventSetup.getData(pTable_);

  HepMC3::GenEvent hepmc_event;
  hepmc_event.set_event_number(event.id().event());
  hepmc_event.set_units(HepMC3::Units::GEV, HepMC3::Units::MM);

  // ---- Event attributes, weights and cross-section (mirror the standard converter) ----
  if (genEventInfoHandle.isValid()) {
    hepmc_event.add_attribute("signal_process_id",
                              std::make_shared<HepMC3::IntAttribute>(genEventInfoHandle->signalProcessID()));
    hepmc_event.add_attribute("event_scale", std::make_shared<HepMC3::DoubleAttribute>(genEventInfoHandle->qScale()));
    hepmc_event.add_attribute("alphaQCD", std::make_shared<HepMC3::DoubleAttribute>(genEventInfoHandle->alphaQCD()));
    hepmc_event.add_attribute("alphaQED", std::make_shared<HepMC3::DoubleAttribute>(genEventInfoHandle->alphaQED()));

    hepmc_event.weights() = genEventInfoHandle->weights();
    if (hepmc_event.weights().empty()) {
      hepmc_event.weights().push_back(1.);
    }

    // PDF info
    const gen::PdfInfo* pdf = genEventInfoHandle->pdf();
    if (pdf != nullptr) {
      const int pdf_id1 = pdf->id.first, pdf_id2 = pdf->id.second;
      const double pdf_x1 = pdf->x.first, pdf_x2 = pdf->x.second;
      const double pdf_scalePDF = pdf->scalePDF;
      const double pdf_xPDF1 = pdf->xPDF.first, pdf_xPDF2 = pdf->xPDF.second;
      HepMC3::GenPdfInfoPtr hepmc_pdfInfo = make_shared<HepMC3::GenPdfInfo>();
      hepmc_pdfInfo->set(pdf_id1, pdf_id2, pdf_x1, pdf_x2, pdf_scalePDF, pdf_xPDF1, pdf_xPDF2);
      hepmc_event.set_pdf_info(hepmc_pdfInfo);
    }
  } else {
    // A graph produced without access to the generator's run / event info still
    // produces a valid HepMC3 event: a unit weight and a dummy cross-section are
    // the convention used by the standard converter when the products are absent.
    hepmc_event.weights().push_back(1.);
  }

  // Cross section, sized to the number of weights
  if (xsec_) {
    if (xsec_->xsecs().size() < hepmc_event.weights().size()) {
      xsec_->set_cross_section(std::vector<double>(hepmc_event.weights().size(), xsec_->xsec(0)),
                              std::vector<double>(hepmc_event.weights().size(), xsec_->xsec_err(0)));
    }
    hepmc_event.set_cross_section(xsec_);
  }

  // ---- Create one HepMC3 GenParticle per kept logical particle ----
  std::unordered_map<uint32_t, HepMC3::GenParticlePtr> particleById;
  particleById.reserve(graph.nParticles());

  HepMC3::GenParticlePtr hepmc_beam1, hepmc_beam2;

  for (uint32_t particleId = 0; particleId < graph.nParticles(); ++particleId) {
    const auto& pdata = graph.particles()[particleId];

    if (!keepParticle(pdata)) {
      continue;
    }

    HepMC3::GenParticlePtr hepmc_particle = std::make_shared<HepMC3::GenParticle>(
        fourVector(pdata.momentum, /*isPosition=*/false), pdata.pdgId, pdata.status);

    // Generated mass: the ParticleDataTable is the same source the standard converter
    // uses, so the truth-graph and the standard HepMC3 records carry the same mass for
    // a given PDG id. Fall back to the stored four-momentum's mass when the particle
    // is not in the table (a quark, an unknown id).
    double particleMass;
    if (pTableData.particle(pdata.pdgId)) {
      particleMass = pTableData.particle(pdata.pdgId)->mass();
    } else {
      particleMass = pdata.momentum.M();
    }
    hepmc_particle->set_generated_mass(particleMass);

    particleById.emplace(particleId, hepmc_particle);

    // While scanning, identify the two incoming beam protons (status 4, forward and
    // pz close to half the cms energy). Used to set the HepMC3 beam particles, the
    // one piece of metadata Rivet strictly requires from the event.
    if (std::abs(pdata.pdgId) == 2212 && pdata.status == 4) {
      if (pdata.momentum.pz() > 0. && !hepmc_beam1) {
        hepmc_beam1 = hepmc_particle;
      } else if (pdata.momentum.pz() < 0. && !hepmc_beam2) {
        hepmc_beam2 = hepmc_particle;
      }
    }
  }

  // ---- Create one HepMC3 GenVertex per kept logical vertex, and wire the topology ----
  for (uint32_t vertexId = 0; vertexId < graph.nVertices(); ++vertexId) {
    const auto& vdata = graph.vertices()[vertexId];

    if (!keepVertex(vdata)) {
      continue;
    }

    HepMC3::GenVertexPtr hepmc_vertex =
        std::make_shared<HepMC3::GenVertex>(fourVector(vdata.position, /*isPosition=*/true));
    hepmc_event.add_vertex(hepmc_vertex);

    for (uint32_t inId : graph.incomingParticles(vertexId)) {
      auto it = particleById.find(inId);
      if (it == particleById.end()) {
        continue;  // particle was skipped
      }
      hepmc_vertex->add_particle_in(it->second);
    }
    for (uint32_t outId : graph.outgoingParticles(vertexId)) {
      auto it = particleById.find(outId);
      if (it == particleById.end()) {
        continue;
      }
      hepmc_vertex->add_particle_out(it->second);
    }
  }

  // ---- Beam particles ----
  // First, if the scan did not identify two status-4 beam protons, look more
  // generally at particles that are particle_in of some vertex but particle_out of
  // no vertex: a graph produced from a HepMC record always has the beam pair here,
  // even when the producer did not stamp status 4 (a legacy HepMC2 record without
  // status 4, or a graph pruned by the selection). The |eta|>5, |pz|>1000 GeV cut
  // used by the standard converter identifies the two pp pair.
  if (!hepmc_beam1 || !hepmc_beam2) {
    for (uint32_t particleId = 0; particleId < graph.nParticles(); ++particleId) {
      const auto& pdata = graph.particles()[particleId];
      auto it = particleById.find(particleId);
      if (it == particleById.end()) {
        continue;
      }
      // Only a real (non-synthetic) particle that has no production vertex can be a
      // beam particle: it must be a graph root with a decay vertex.
      if (pdata.isSynthetic()) {
        continue;
      }
      if (!graph.productionVertices(particleId).empty()) {
        continue;
      }
      if (graph.decayVertices(particleId).empty()) {
        continue;
      }
      const auto& p = pdata.momentum;
      if (std::abs(p.eta()) <= 5. || std::abs(p.pz()) <= 1000.) {
        continue;
      }
      if (p.pz() > 0. && !hepmc_beam1) {
        hepmc_beam1 = it->second;
      } else if (p.pz() < 0. && !hepmc_beam2) {
        hepmc_beam2 = it->second;
      }
      if (hepmc_beam1 && hepmc_beam2) {
        break;
      }
    }
  }

  // If the graph carried the beam protons, the HepMC3 record already has the
  // interaction vertex with them as particle_in, faithfully reproducing the
  // generator record. This is the case the truth graph is built for.
  if (hepmc_beam1 && hepmc_beam2) {
    hepmc_event.set_beam_particles(hepmc_beam1, hepmc_beam2);
  } else if (!particleById.empty()) {
    // No beam pair was found: the input graph either has no generator record at
    // all (a particle-gun sample, or a graph pruned by a selection preset). Synthesize
    // a minimal interaction vertex with two dummy beam protons, so the HepMC3 record
    // the standard converter produces can still be reproduced, with the incident
    // protons reconstructed from the configured cmEnergy.
    //
    // The motherless outgoing particles are attached to one of the two incident
    // vertices (pz > 0 -> +z, pz < 0 -> -z), matching the rule the standard converter
    // uses to attach a record without incoming partons to a beam remnant side.
    const double beamEnergy = cmEnergy_ / 2.;
    const HepMC3::FourVector nullVtx(0., 0., 0., 0.);
    HepMC3::FourVector beamPlus(0., 0., +beamEnergy, beamEnergy);
    HepMC3::FourVector beamMinus(0., 0., -beamEnergy, beamEnergy);

    HepMC3::GenParticlePtr dummyBeam1 = std::make_shared<HepMC3::GenParticle>(beamPlus, 2212, 4);
    HepMC3::GenParticlePtr dummyBeam2 = std::make_shared<HepMC3::GenParticle>(beamMinus, 2212, 4);

    HepMC3::GenVertexPtr vtx1 = std::make_shared<HepMC3::GenVertex>(nullVtx);
    HepMC3::GenVertexPtr vtx2 = std::make_shared<HepMC3::GenVertex>(nullVtx);
    hepmc_event.add_vertex(vtx1);
    hepmc_event.add_vertex(vtx2);
    vtx1->add_particle_in(dummyBeam1);
    vtx2->add_particle_in(dummyBeam2);

    for (uint32_t particleId = 0; particleId < graph.nParticles(); ++particleId) {
      auto it = particleById.find(particleId);
      if (it == particleById.end()) {
        continue;
      }
      // Skip particles that already have a production vertex: their topology has
      // already been wired and they belong to their real vertex, not the dummy one.
      if (!graph.productionVertices(particleId).empty()) {
        continue;
      }
      // A particle that's already particle_in somewhere also has its topology
      // attached already: it's a beam particle of the original record. Skip it.
      bool alreadyInSomeVertex = false;
      for (uint32_t decayVtx : graph.decayVertices(particleId)) {
        for (uint32_t incomingId : graph.incomingParticles(decayVtx)) {
          if (incomingId == particleId) {
            alreadyInSomeVertex = true;
            break;
          }
        }
        if (alreadyInSomeVertex) {
          break;
      }
      }
      if (alreadyInSomeVertex) {
        continue;
      }

      if (graph.particles()[particleId].momentum.pz() > 0.) {
        vtx1->add_particle_out(it->second);
      } else {
        vtx2->add_particle_out(it->second);
      }
    }

    hepmc_event.set_beam_particles(dummyBeam1, dummyBeam2);
  }

  // ---- Write the HepMC3 record to the text file, if requested ----
  if (writeHepMC_) {
    writer_->write_event(hepmc_event);
  }

  if (verbosity_ > 0) {
    edm::LogVerbatim("TruthGraph2HepMCConverter")
        << "truth::Graph -> HepMC3: " << particleById.size() << " / " << graph.nParticles()
        << " particles emitted, " << hepmc_event.vertices().size() << " vertices, "
        << (hepmc_event.beams().size() == 2 ? "beam pair set"
                                             : "no beam pair (" + std::to_string(hepmc_event.beams().size()) + ")")
        << ", event " << event.id().event();
  }

  // ---- Finalize and put the product into the event ----
  auto hepmc_product = std::make_unique<edm::HepMC3Product>(hepmc_event);
  event.put(std::move(hepmc_product), "unsmeared");
}

DEFINE_FWK_MODULE(TruthGraph2HepMCConverter);
