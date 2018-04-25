// -*- C++ -*-
//
// Package:    PhysicsTools/NanoAOD
// Class:      IsoValueMapProducer
//
/**\class IsoValueMapProducer IsoValueMapProducer.cc PhysicsTools/NanoAOD/plugins/IsoValueMapProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco Peruzzi
//         Created:  Mon, 04 Sep 2017 22:43:53 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"

//
// class declaration
//

template <typename T>
class IsoValueMapProducer : public edm::global::EDProducer<> {
   public:
  explicit IsoValueMapProducer(const edm::ParameterSet &iConfig):
    src_(consumes<edm::View<T>>(iConfig.getParameter<edm::InputTag>("src"))),
    relative_(iConfig.getParameter<bool>("relative"))
  {
    if ((typeid(T) == typeid(pat::Muon)) || (typeid(T) == typeid(pat::Electron)) || typeid(T) == typeid(pat::IsolatedTrack)) {
      produces<edm::ValueMap<float>>("miniIsoChg");
      produces<edm::ValueMap<float>>("miniIsoAll");
      ea_miniiso_.reset(new EffectiveAreas((iConfig.getParameter<edm::FileInPath>("effectiveAreas")).fullPath()));
      rho_miniiso_ = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));
      vertexes_ = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));
      pfCandidates_ = consumes<std::vector<reco::PFCandidate>>(iConfig.getParameter<edm::InputTag>("pfCandidates"));
    }
    if ((typeid(T) == typeid(pat::Electron)) or (typeid(T) == typeid(pat::Muon))){
      ea_pfiso_.reset(new EffectiveAreas((iConfig.getParameter<edm::FileInPath>("effectiveAreas")).fullPath()));
      rho_pfiso_ = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));
      produces<edm::ValueMap<float>>("PFIsoChg");
      produces<edm::ValueMap<float>>("PFIsoAll");
    }
    else if ((typeid(T) == typeid(pat::Photon))) {
      produces<edm::ValueMap<float>>("PFIsoChg");
      produces<edm::ValueMap<float>>("PFIsoAll");
      mapIsoChg_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mapIsoChg"));
      mapIsoNeu_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mapIsoNeu"));
      mapIsoPho_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mapIsoPho"));
      ea_pfiso_chg_.reset(new EffectiveAreas((iConfig.getParameter<edm::FileInPath>("EAFile_PFIso_Chg")).fullPath()));
      ea_pfiso_neu_.reset(new EffectiveAreas((iConfig.getParameter<edm::FileInPath>("EAFile_PFIso_Neu")).fullPath()));
      ea_pfiso_pho_.reset(new EffectiveAreas((iConfig.getParameter<edm::FileInPath>("EAFile_PFIso_Pho")).fullPath()));
      rho_pfiso_ = consumes<double>(iConfig.getParameter<edm::InputTag>("rho_PFIso"));
    }
  }
  ~IsoValueMapProducer() override {}

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

      // ----------member data ---------------------------

  edm::EDGetTokenT<edm::View<T>> src_;
  bool relative_;
  edm::EDGetTokenT<double> rho_miniiso_;
  edm::EDGetTokenT<reco::VertexCollection> vertexes_;
  edm::EDGetTokenT<double> rho_pfiso_;
  edm::EDGetTokenT<edm::ValueMap<float>> mapIsoChg_;
  edm::EDGetTokenT<edm::ValueMap<float>> mapIsoNeu_;
  edm::EDGetTokenT<edm::ValueMap<float>> mapIsoPho_;
  edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfCandidates_;
  std::unique_ptr<EffectiveAreas> ea_miniiso_;
  std::unique_ptr<EffectiveAreas> ea_pfiso_;
  std::unique_ptr<EffectiveAreas> ea_pfiso_chg_;
  std::unique_ptr<EffectiveAreas> ea_pfiso_neu_;
  std::unique_ptr<EffectiveAreas> ea_pfiso_pho_;
  float getEtaForEA(const T*) const;
  void doMiniIso(edm::Event&) const;
  void doPFIsoEle(edm::Event&) const;
  void doPFIsoMuon(edm::Event&) const;
  void doPFIsoPho(edm::Event&) const;

  double chargedSum(const T&, edm::Handle<reco::PFCandidateCollection>, double, edm::Handle<reco::VertexCollection> PVs) const;
  double neutralSum(const T&, edm::Handle<reco::PFCandidateCollection>, double) const;
  double photonSum(const T&, edm::Handle<reco::PFCandidateCollection>, double) const;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

template<typename T> float IsoValueMapProducer<T>::getEtaForEA(const T *obj) const{
  return obj->eta();
}
template<> float IsoValueMapProducer<pat::Electron>::getEtaForEA(const pat::Electron *el) const{
  return el->superCluster()->eta();
}
template<> float IsoValueMapProducer<pat::Photon>::getEtaForEA(const pat::Photon *ph) const{
  return ph->superCluster()->eta();
}

template <typename T>
void
IsoValueMapProducer<T>::produce(edm::StreamID streamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const
{

  if ((typeid(T) == typeid(pat::Muon)) || (typeid(T) == typeid(pat::Electron)) || typeid(T) == typeid(pat::IsolatedTrack)) { doMiniIso(iEvent); };
  if ((typeid(T) == typeid(pat::Electron))) { doPFIsoEle(iEvent); }
  if ((typeid(T) == typeid(pat::Muon))) { doPFIsoMuon(iEvent); }
  if ((typeid(T) == typeid(pat::Photon))) { doPFIsoPho(iEvent); }

}


template<typename T> double IsoValueMapProducer<T>::chargedSum(const T& muon, edm::Handle<reco::PFCandidateCollection> pfcands, double coneSize, edm::Handle<reco::VertexCollection> PVs) const{
  double deadcone = 0.0001;

  double iso_ch = 0;
  for(const reco::PFCandidate &pfc : *pfcands){
    if(fabs(pfc.pdgId()) < 7)  continue;
    if(pfc.charge()==0)        continue;


    auto vtxLead  = PVs->begin();
    auto vtxClose = PVs->begin();
    for(auto vtx = PVs->begin(); vtx != PVs->end(); ++vtx){
      if(fabs(pfc.trackRef()->dz(vtx->position())) < fabs(pfc.trackRef()->dz(vtxClose->position()))) vtxClose = vtx;
    }
    if(vtxClose != vtxLead)    continue; // AOD equivalent of fromPV > 1

    if(fabs(pfc.pdgId())!=211) continue;

    double dr = deltaR(pfc, muon);
    if(dr < coneSize and dr > deadcone) iso_ch += pfc.pt();
  }
  return iso_ch;
}

template<typename T> double IsoValueMapProducer<T>::neutralSum(const T& muon, edm::Handle<reco::PFCandidateCollection> pfcands, double coneSize) const{
  double deadcone = 0.01;

  double iso_nh = 0;
  for(const reco::PFCandidate &pfc : *pfcands){
    if(fabs(pfc.pdgId()) < 7)  continue;
    if(pfc.charge()!=0)        continue;
    if(fabs(pfc.pdgId())!=130) continue;
    if(pfc.pt()<0.5)           continue;

    double dr = deltaR(pfc, muon);
    if(dr < coneSize and dr > deadcone) iso_nh += pfc.pt();
  }
  return iso_nh;
}

template<typename T> double IsoValueMapProducer<T>::photonSum(const T& muon, edm::Handle<reco::PFCandidateCollection> pfcands, double coneSize) const{
  double deadcone = 0.01; 

  double iso_ph = 0;
  for(const reco::PFCandidate &pfc : *pfcands){
    if(fabs(pfc.pdgId()) < 7)  continue;
    if(pfc.charge()!=0)        continue;
    if(fabs(pfc.pdgId())!=22)  continue;
    if(pfc.pt()<0.5)           continue;

    double dr = deltaR(pfc, muon);
    if(dr < coneSize and dr > deadcone) iso_ph += pfc.pt();
  }
  return iso_ph;
}


template<typename T>
void
IsoValueMapProducer<T>::doMiniIso(edm::Event& iEvent) const{

  edm::Handle<edm::View<T>> src;
  iEvent.getByToken(src_, src);
  edm::Handle<double> rho;
  iEvent.getByToken(rho_miniiso_,rho);
  edm::Handle<std::vector<reco::PFCandidate>> pfCandidates;
  iEvent.getByToken(pfCandidates_, pfCandidates);

  edm::Handle<reco::VertexCollection> PVs;
  iEvent.getByToken(vertexes_, PVs);


  unsigned int nInput = src->size();

  std::vector<float> miniIsoChg, miniIsoAll;
  miniIsoChg.reserve(nInput);
  miniIsoAll.reserve(nInput);

  for (const auto & obj : *src) {
    auto iso = obj.miniPFIsolation();
    auto ea = ea_miniiso_->getEffectiveArea(fabs(getEtaForEA(&obj)));
    float R = 10.0/std::min(std::max(obj.pt(), 50.0),200.0);
    auto chg = typeid(T) != typeid(pat::Muon) ? iso.chargedHadronIso() : chargedSum(obj, pfCandidates, R, PVs);
    auto neu = typeid(T) != typeid(pat::Muon) ? iso.neutralHadronIso() : neutralSum(obj, pfCandidates, R);
    auto pho = typeid(T) != typeid(pat::Muon) ? iso.photonIso()        : photonSum(obj, pfCandidates, R);
    ea *= std::pow(R/0.3,2);
    float scale = relative_ ? 1.0/obj.pt() : 1;
    miniIsoChg.push_back(scale*chg);
    miniIsoAll.push_back(scale*(chg+std::max(0.0,neu+pho-(*rho)*ea)));
  }

  std::unique_ptr<edm::ValueMap<float>> miniIsoChgV(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerChg(*miniIsoChgV);
  fillerChg.insert(src,miniIsoChg.begin(),miniIsoChg.end());
  fillerChg.fill();
  std::unique_ptr<edm::ValueMap<float>> miniIsoAllV(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerAll(*miniIsoAllV);
  fillerAll.insert(src,miniIsoAll.begin(),miniIsoAll.end());
  fillerAll.fill();

  iEvent.put(std::move(miniIsoChgV),"miniIsoChg");
  iEvent.put(std::move(miniIsoAllV),"miniIsoAll");
}

template<>
void
IsoValueMapProducer<pat::Photon>::doMiniIso(edm::Event& iEvent) const {}


template<typename T> void IsoValueMapProducer<T>::doPFIsoEle(edm::Event& iEvent) const {}
template<typename T> void IsoValueMapProducer<T>::doPFIsoMuon(edm::Event& iEvent) const {}

template<> void IsoValueMapProducer<pat::Electron>::doPFIsoEle(edm::Event& iEvent) const{
  edm::Handle<edm::View<pat::Electron>> src;
  iEvent.getByToken(src_, src);
  edm::Handle<double> rho;
  iEvent.getByToken(rho_pfiso_,rho);
  edm::Handle<std::vector<reco::PFCandidate>> pfCandidates;
  iEvent.getByToken(pfCandidates_, pfCandidates);

  unsigned int nInput = src->size();

  std::vector<float> PFIsoChg, PFIsoAll;
  PFIsoChg.reserve(nInput);
  PFIsoAll.reserve(nInput);

  for (const auto & obj : *src) {
    auto iso = obj.pfIsolationVariables();
    auto chg = iso.sumChargedHadronPt;
    auto neu = iso.sumNeutralHadronEt;
    auto pho = iso.sumPhotonEt;
    auto ea = ea_pfiso_->getEffectiveArea(fabs(getEtaForEA(&obj)));
    float scale = relative_ ? 1.0/obj.pt() : 1;
    PFIsoChg.push_back(scale*chg);
    PFIsoAll.push_back(scale*(chg+std::max(0.0,neu+pho-(*rho)*ea)));
  }

  std::unique_ptr<edm::ValueMap<float>> PFIsoChgV(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerChg(*PFIsoChgV);
  fillerChg.insert(src,PFIsoChg.begin(),PFIsoChg.end());
  fillerChg.fill();
  std::unique_ptr<edm::ValueMap<float>> PFIsoAllV(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerAll(*PFIsoAllV);
  fillerAll.insert(src,PFIsoAll.begin(),PFIsoAll.end());
  fillerAll.fill();

  iEvent.put(std::move(PFIsoChgV),"PFIsoChg");
  iEvent.put(std::move(PFIsoAllV),"PFIsoAll");
}

template<> void IsoValueMapProducer<pat::Muon>::doPFIsoMuon(edm::Event& iEvent) const{
  edm::Handle<edm::View<pat::Muon>> src;
  iEvent.getByToken(src_, src);
  edm::Handle<double> rho;
  iEvent.getByToken(rho_pfiso_,rho);
  edm::Handle<std::vector<reco::PFCandidate>> pfCandidates;
  iEvent.getByToken(pfCandidates_, pfCandidates);

  unsigned int nInput = src->size();

  std::vector<float> PFIsoChg, PFIsoAll;
  PFIsoChg.reserve(nInput);
  PFIsoAll.reserve(nInput);

  for (const auto & obj : *src) {
    auto iso = obj.pfIsolationR03();
    auto chg = iso.sumChargedHadronPt;
    auto neu = iso.sumNeutralHadronEt;
    auto pho = iso.sumPhotonEt;
    auto ea = ea_pfiso_->getEffectiveArea(fabs(obj.eta()));
    float scale = relative_ ? 1.0/obj.pt() : 1;
    PFIsoChg.push_back(scale*chg);
    PFIsoAll.push_back(scale*(chg+std::max(0.0,neu+pho-(*rho)*ea)));
  }

  std::unique_ptr<edm::ValueMap<float>> PFIsoChgV(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerChg(*PFIsoChgV);
  fillerChg.insert(src,PFIsoChg.begin(),PFIsoChg.end());
  fillerChg.fill();
  std::unique_ptr<edm::ValueMap<float>> PFIsoAllV(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerAll(*PFIsoAllV);
  fillerAll.insert(src,PFIsoAll.begin(),PFIsoAll.end());
  fillerAll.fill();

  iEvent.put(std::move(PFIsoChgV),"PFIsoChg");
  iEvent.put(std::move(PFIsoAllV),"PFIsoAll");
}

template<typename T>
void
IsoValueMapProducer<T>::doPFIsoPho(edm::Event& iEvent) const {}

template<>
void
IsoValueMapProducer<pat::Photon>::doPFIsoPho(edm::Event& iEvent) const {

  edm::Handle<edm::View<pat::Photon>> src;
  iEvent.getByToken(src_, src);
  edm::Handle<double> rho;
  iEvent.getByToken(rho_pfiso_,rho);
  edm::Handle<edm::ValueMap<float> > mapIsoChg;
  iEvent.getByToken(mapIsoChg_, mapIsoChg);
  edm::Handle<edm::ValueMap<float> > mapIsoNeu;
  iEvent.getByToken(mapIsoNeu_, mapIsoNeu);
  edm::Handle<edm::ValueMap<float> > mapIsoPho;
  iEvent.getByToken(mapIsoPho_, mapIsoPho);

  unsigned int nInput = src->size();

  std::vector<float> PFIsoChg, PFIsoAll;
  PFIsoChg.reserve(nInput);
  PFIsoAll.reserve(nInput);

  for (unsigned int i=0; i<nInput; i++){
    auto obj = src->ptrAt(i);
    auto chg = (*mapIsoChg)[obj];
    auto neu = (*mapIsoNeu)[obj];
    auto pho = (*mapIsoPho)[obj];
    auto ea_chg = ea_pfiso_chg_->getEffectiveArea(fabs(getEtaForEA(obj.get())));
    auto ea_neu = ea_pfiso_neu_->getEffectiveArea(fabs(getEtaForEA(obj.get())));
    auto ea_pho = ea_pfiso_pho_->getEffectiveArea(fabs(getEtaForEA(obj.get())));
    float scale = relative_ ? 1.0/obj->pt() : 1;
    PFIsoChg.push_back(scale*std::max(0.0,chg-(*rho)*ea_chg));
    PFIsoAll.push_back(PFIsoChg.back()+scale*(std::max(0.0,neu-(*rho)*ea_neu)+std::max(0.0,pho-(*rho)*ea_pho)));
  }

  std::unique_ptr<edm::ValueMap<float>> PFIsoChgV(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerChg(*PFIsoChgV);
  fillerChg.insert(src,PFIsoChg.begin(),PFIsoChg.end());
  fillerChg.fill();
  std::unique_ptr<edm::ValueMap<float>> PFIsoAllV(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerAll(*PFIsoAllV);
  fillerAll.insert(src,PFIsoAll.begin(),PFIsoAll.end());
  fillerAll.fill();

  iEvent.put(std::move(PFIsoChgV),"PFIsoChg");
  iEvent.put(std::move(PFIsoAllV),"PFIsoAll");

}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template <typename T>
void
IsoValueMapProducer<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src")->setComment("input physics object collection");
  desc.add<bool>("relative")->setComment("compute relative isolation instead of absolute one");
  if ((typeid(T) == typeid(pat::Muon)) || (typeid(T) == typeid(pat::Electron)) || typeid(T) == typeid(pat::IsolatedTrack)) {
    desc.add<edm::FileInPath>("effectiveAreas")->setComment("txt file containing effective areas to be used for mini-isolation pileup subtraction");
    desc.add<edm::InputTag>("rho")->setComment("rho to be used for effective-area based mini-isolation pileup subtraction");
    desc.add<edm::InputTag>("pfCandidates")->setComment("pfCandidates to calculate isolation sums for mini-isolation");
  }
  if ((typeid(T) == typeid(pat::Electron))) {
    desc.add<edm::FileInPath>("EAFile_PFIso")->setComment("txt file containing effective areas to be used for PF-isolation pileup subtraction for electrons");
    desc.add<edm::InputTag>("rho_PFIso")->setComment("rho to be used for effective-area based PF-isolation pileup subtraction for electrons");
  }
  if ((typeid(T) == typeid(pat::Photon))) {
    desc.add<edm::InputTag>("mapIsoChg")->setComment("input charged PF isolation calculated in VID for photons");
    desc.add<edm::InputTag>("mapIsoNeu")->setComment("input neutral PF isolation calculated in VID for photons");
    desc.add<edm::InputTag>("mapIsoPho")->setComment("input photon PF isolation calculated in VID for photons");
    desc.add<edm::FileInPath>("EAFile_PFIso_Chg")->setComment("txt file containing effective areas to be used for charged PF-isolation pileup subtraction for photons");
    desc.add<edm::FileInPath>("EAFile_PFIso_Neu")->setComment("txt file containing effective areas to be used for neutral PF-isolation pileup subtraction for photons");
    desc.add<edm::FileInPath>("EAFile_PFIso_Pho")->setComment("txt file containing effective areas to be used for photon PF-isolation pileup subtraction for photons");
    desc.add<edm::InputTag>("rho_PFIso")->setComment("rho to be used for effective-area based PF-isolation pileup subtraction for photons");
  }
  std::string modname;
  if (typeid(T) == typeid(pat::Muon)) modname+="Muon";
  else if (typeid(T) == typeid(pat::Electron)) modname+="Ele";
  else if (typeid(T) == typeid(pat::Photon)) modname+="Pho";
  else if (typeid(T) == typeid(pat::IsolatedTrack)) modname+="IsoTrack";
  modname+="IsoValueMapProducer";
  descriptions.add(modname,desc);
}


typedef IsoValueMapProducer<pat::Muon> MuonIsoValueMapProducer;
typedef IsoValueMapProducer<pat::Electron> EleIsoValueMapProducer;
typedef IsoValueMapProducer<pat::Photon> PhoIsoValueMapProducer;
typedef IsoValueMapProducer<pat::IsolatedTrack> IsoTrackIsoValueMapProducer;

//define this as a plug-in
DEFINE_FWK_MODULE(MuonIsoValueMapProducer);
DEFINE_FWK_MODULE(EleIsoValueMapProducer);
DEFINE_FWK_MODULE(PhoIsoValueMapProducer);
DEFINE_FWK_MODULE(IsoTrackIsoValueMapProducer);

