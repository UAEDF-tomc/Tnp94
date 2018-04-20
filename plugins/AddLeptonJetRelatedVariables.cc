// system include files
#include <memory>
#include <cmath>
#include <TLorentzVector.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include <DataFormats/MuonReco/interface/Muon.h>
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/Math/interface/deltaR.h"

using namespace std;


//
// class declaration
//

class AddLeptonJetRelatedVariables : public edm::EDProducer{
public:
  explicit AddLeptonJetRelatedVariables(const edm::ParameterSet&);
  ~AddLeptonJetRelatedVariables();

  typedef std::vector< edm::FwdPtr<reco::PFCandidate> > PFCollection;
  
private:
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  

  // ----------auxiliary functions -------------------  
  // ----------member data ---------------------------
  edm::EDGetTokenT<reco::PFJetCollection> jetCollectionTag_;
  edm::EDGetTokenT<reco::CandidateView> leptonCollectionTag_;
  edm::EDGetTokenT<reco::JetCorrector> tagL1Corrector_;
  edm::EDGetTokenT<reco::JetCorrector> tagL1L2L3ResCorrector_;
  edm::EDGetTokenT<reco::VertexCollection> vertexes_;
  edm::EDGetTokenT<reco::JetTagCollection> bTagCollectionTag_;
  edm::EDGetTokenT<reco::JetTagCollection> bTagCollectionTag2_;
  edm::EDGetTokenT<reco::JetTagCollection> bTagCollectionTag3_;


  template<typename Hand, typename T>
  void storeMap(edm::Event &iEvent, const Hand & handle, const std::vector<T> & values, const std::string    & label) const ; 
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
AddLeptonJetRelatedVariables::AddLeptonJetRelatedVariables(const edm::ParameterSet& iConfig){
  edm::InputTag jetcollection = iConfig.getParameter<edm::InputTag>("RawJetCollection");
  jetCollectionTag_ = consumes<reco::PFJetCollection>(jetcollection);

  edm::InputTag leptoncollection = iConfig.getParameter<edm::InputTag>("LeptonCollection");
  leptonCollectionTag_ = consumes<reco::CandidateView>(leptoncollection);

  edm::InputTag l1Cortag = iConfig.getParameter<edm::InputTag>("L1Corrector");
  tagL1Corrector_ = consumes<reco::JetCorrector>(l1Cortag);

  edm::InputTag l1l2l3ResCortag = iConfig.getParameter<edm::InputTag>("L1L2L3ResCorrector");
  tagL1L2L3ResCorrector_ = consumes<reco::JetCorrector>(l1l2l3ResCortag);

  bTagCollectionTag_  = consumes<reco::JetTagCollection>(edm::InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
  bTagCollectionTag2_ = consumes<reco::JetTagCollection>(edm::InputTag("pfDeepCSVJetTags:probb"));
  bTagCollectionTag3_ = consumes<reco::JetTagCollection>(edm::InputTag("pfDeepCSVJetTags:probbb"));
  
  vertexes_ = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));

  //now do what ever initialization is needed
  produces<edm::ValueMap<float> >("JetPtRatio");
  produces<edm::ValueMap<float> >("JetPtRel");
  produces<edm::ValueMap<float> >("JetNDauCharged");
  produces<edm::ValueMap<float> >("JetBTagCSV");
  produces<edm::ValueMap<float> >("JetDeepBTagCSV");
}

AddLeptonJetRelatedVariables::~AddLeptonJetRelatedVariables()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//




// ------------ method called for each event  ------------
void 
AddLeptonJetRelatedVariables::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
  using namespace edm;

  std::cout << std::endl << "Run: " << iEvent.id().run() << "\tLumi: " << iEvent.id().luminosityBlock() << "\tEvent: " << iEvent.id().event() << std::endl;
  edm::Handle<reco::PFJetCollection> jets;       
  iEvent.getByToken (jetCollectionTag_, jets);    

  edm::Handle<reco::CandidateView> leptons;               
  iEvent.getByToken (leptonCollectionTag_, leptons);       

  edm::Handle<reco::JetCorrector> correctorL1L2L3Res;
  iEvent.getByToken(tagL1L2L3ResCorrector_, correctorL1L2L3Res);

  edm::Handle<reco::JetCorrector> correctorL1;
  iEvent.getByToken(tagL1Corrector_, correctorL1);
  
  edm::Handle<reco::JetTagCollection> bTagHandle;  iEvent.getByToken(bTagCollectionTag_,  bTagHandle);  const reco::JetTagCollection & bTags  = *(bTagHandle.product());
  edm::Handle<reco::JetTagCollection> bTagHandle2; iEvent.getByToken(bTagCollectionTag2_, bTagHandle2); const reco::JetTagCollection & bTags2 = *(bTagHandle2.product());
  edm::Handle<reco::JetTagCollection> bTagHandle3; iEvent.getByToken(bTagCollectionTag3_, bTagHandle3); const reco::JetTagCollection & bTags3 = *(bTagHandle3.product());

  edm::Handle<reco::VertexCollection> PVs;
  iEvent.getByToken(vertexes_, PVs);
  reco::VertexRef PV(PVs.id());
  reco::VertexRefProd PVRefProd(PVs);

  //Output
  std::vector<float> ptratio,ptrel,nchargeddaughers,btagcsv,deepbtagcsv;
  for(auto icand = leptons->begin(); icand != leptons->end(); ++ icand){

    // Find closest selected jet
    std::vector<pat::Jet> selectedJetsAll;
    std::vector<float> csvs, deepcsvs;
    for(auto jet = jets->begin(); jet != jets->end(); ++jet){
      if(jet->pt() > 5 && fabs( jet->eta() ) < 3) selectedJetsAll.push_back(*jet);
    }
    unsigned closestIndex = 0;
    for(unsigned j = 1; j < selectedJetsAll.size(); ++j){
      if(reco::deltaR(selectedJetsAll[j], *icand) < reco::deltaR(selectedJetsAll[closestIndex], *icand)) closestIndex = j;
    }

    const pat::Jet& jet = selectedJetsAll[closestIndex];
    if(selectedJetsAll.size() == 0 || reco::deltaR(jet, *icand) > 0.4){
      ptratio.push_back(1);
      ptrel.push_back(0);
      btagcsv.push_back(0);
      deepbtagcsv.push_back(0);
      nchargeddaughers.push_back(0);
    } else {

      float bjet(0);
      float bjet2(0);
      float bjet3(0);
      int nchdaugthers(0);    
      
      for (reco::JetTagCollection::const_iterator tagI = bTags.begin(); tagI != bTags.end(); ++tagI) {
        if(fabs(tagI->first->pt() - jet.pt()) < 0.001) bjet = tagI->second;
      }

      for (reco::JetTagCollection::const_iterator tagI = bTags2.begin(); tagI != bTags2.end(); ++tagI) {
        if(fabs(tagI->first->pt() - jet.pt()) < 0.001) bjet2 = tagI->second;
      }

      for (reco::JetTagCollection::const_iterator tagI = bTags3.begin(); tagI != bTags3.end(); ++tagI) {
        if(fabs(tagI->first->pt() - jet.pt()) < 0.001) bjet3 = tagI->second;
      }

      btagcsv.push_back(bjet);
      deepbtagcsv.push_back(bjet2+bjet3);

      PV = reco::VertexRef(PVs, 0);
      math::XYZPoint PVpos = PV->position();
      for(auto part : jet.getJetConstituentsQuick()){
          auto daughter = static_cast<const reco::PFCandidate*>(part);
          reco::TrackRef daughterTrack = daughter->trackRef();
          if(not daughterTrack.isNonnull())      continue;
          if(reco::deltaR(jet, *daughter) > 0.4) continue;

          bool goodTrack = daughterTrack->pt() > 1 and
                           daughterTrack->charge() != 0 and
                           daughterTrack->hitPattern().numberOfValidHits() > 7 and
                           daughterTrack->hitPattern().numberOfValidPixelHits() > 1 and
                           daughterTrack->normalizedChi2() < 5 and
                           fabs(daughterTrack->dz(PVpos)) < 17 and
                           fabs(daughterTrack->dxy(PVpos)) < 0.2;
          if(not goodTrack) continue;

          auto vtxLead  = PVs->begin();
          auto vtxClose = PVs->begin();
          for(auto vtx = PVs->begin(); vtx != PVs->end(); ++vtx){
            if(fabs(daughterTrack->dz(vtx->position())) < fabs(daughterTrack->dz(vtxClose->position()))) vtxClose = vtx;
          }
          if(vtxClose == vtxLead) ++nchdaugthers; // AOD equivalent of fromPV > 1
      }

      nchargeddaughers.push_back(nchdaugthers);

      float jecL1L2L3Res = correctorL1L2L3Res->correction(jet);
      float jecL1       = correctorL1->correction(jet);

      float JEC         = jecL1L2L3Res/jecL1;
      auto  l1Jet       = jet.p4()*jecL1;
      auto  l           = icand->p4();
      auto  lepAwareJet = (l1Jet - l)*JEC + l;

      TLorentzVector lV(l.Px(), l.Py(), l.Pz(), l.E());
      TLorentzVector jV(lepAwareJet.Px(), lepAwareJet.Py(), lepAwareJet.Pz(), lepAwareJet.E());
      ptratio.push_back(l.Pt()/lepAwareJet.Pt());
      ptrel.push_back(lV.Perp((jV - lV).Vect()));

      if(true and icand->pt() > 10){ // For debugging
        std::cout << "  **** Muon (pt/eta,phi): " << icand->pt() << ", " << icand->eta() << "\t" << icand->phi() << std::endl;
        std::cout << "        jet (pt/eta,phi): " << lepAwareJet.pt() << ", " << jet.eta() << "\t" << jet.phi() << std::endl;
        std::cout << "                 deepcsv: " << bjet2+bjet3 << "(" << bjet2 << "+" << bjet3 << ")" << std::endl;
        std::cout << "                     csv: " << bjet << std::endl;
        std::cout << "            nchdaugthers: " << nchdaugthers << std::endl;
        std::cout << "                   ptrel: " << ptrel.back() << std::endl;
        std::cout << "                 ptratio: " << ptratio.back() << std::endl;
      }
    }

    //    printf ("muon: %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\n", mu.pt(),mu.eta(),mu.phi(), ptrel.back(), ptratio.back(), nchargeddaughers.back(), btagcsv.back());


  }//end muon loop
  
  
  /// Filling variables previously computed
  storeMap(iEvent, leptons, ptratio, "JetPtRatio");
  storeMap(iEvent, leptons, ptrel, "JetPtRel");
  storeMap(iEvent, leptons, nchargeddaughers, "JetNDauCharged");
  storeMap(iEvent, leptons, btagcsv, "JetBTagCSV");
  storeMap(iEvent, leptons, deepbtagcsv, "JetDeepBTagCSV");

}



// ------------ method called once each job just before starting event loop  ------------
void 
AddLeptonJetRelatedVariables::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AddLeptonJetRelatedVariables::endJob() 
{
}

template<typename Hand, typename T>
void
AddLeptonJetRelatedVariables::storeMap(edm::Event &iEvent,
		    const Hand & handle,
		    const std::vector<T> & values,
		    const std::string    & label) const {
  using namespace edm; using namespace std;
  unique_ptr<ValueMap<T> > valMap(new ValueMap<T>());
  typename edm::ValueMap<T>::Filler filler(*valMap);
  filler.insert(handle, values.begin(), values.end());
  filler.fill();
  iEvent.put(std::move(valMap), label);
}


//define this as a plug-in
DEFINE_FWK_MODULE(AddLeptonJetRelatedVariables);
