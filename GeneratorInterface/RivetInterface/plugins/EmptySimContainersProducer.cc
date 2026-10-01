// Emits empty edm::SimTrackContainer and edm::SimVertexContainer products.
//
// The TruthGraphProducer reads SimTracks/SimVertices via evt.get(), which
// throws if the products are absent from the event. A generator-only job
// (a Pythia8 GEN step with no Geant4 SIM step) does not produce those
// products, so the TruthGraphProducer cannot run there even though it
// handles an empty SimTrack/SimVertex container perfectly well (it just
// builds a GEN-only TruthGraph).
//
// This tiny producer bridges that gap: it puts empty SimTrack/SimVertex
// containers into the event under configurable labels, so the
// TruthGraphProducer's evt.get() finds them and the GEN-only truth graph
// can be built. It has no other purpose.

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

class EmptySimContainersProducer : public edm::stream::EDProducer<> {
public:
  explicit EmptySimContainersProducer(const edm::ParameterSet& pset)
      : simTracksLabel_(pset.getParameter<edm::InputTag>("simTracks")),
        simVerticesLabel_(pset.getParameter<edm::InputTag>("simVertices")) {
    produces<edm::SimTrackContainer>(simTracksLabel_.instance());
    produces<edm::SimVertexContainer>(simVerticesLabel_.instance());
  }

  void produce(edm::Event& evt, const edm::EventSetup&) override {
    evt.put(std::make_unique<edm::SimTrackContainer>(), simTracksLabel_.instance());
    evt.put(std::make_unique<edm::SimVertexContainer>(), simVerticesLabel_.instance());
  }

private:
  edm::InputTag simTracksLabel_;
  edm::InputTag simVerticesLabel_;
};

DEFINE_FWK_MODULE(EmptySimContainersProducer);
