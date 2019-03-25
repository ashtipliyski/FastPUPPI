// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include <algorithm>

class L1MultPuppiWeightProducer : public edm::global::EDProducer<>  {
public:
  explicit L1MultPuppiWeightProducer(const edm::ParameterSet&);
  ~L1MultPuppiWeightProducer();

private:
  virtual void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

  int getMatch(int i, l1t::PFCandidate pfCand, std::vector<l1t::PFCandidate> pupCands) const;

  edm::EDGetTokenT<std::vector<l1t::PFCandidate>> pfToken_;
  edm::EDGetTokenT<std::vector<l1t::PFCandidate>> pupToken_;
  edm::EDGetTokenT<std::vector<l1t::PFCandidate>> pupToken1_;
  edm::EDGetTokenT<std::vector<l1t::PFCandidate>> pupToken2_;
  edm::EDGetTokenT<std::vector<l1t::PFCandidate>> pupToken3_;
  edm::EDGetTokenT<std::vector<l1t::PFCandidate>> pupToken4_;

  unsigned int top_m_;
};

L1MultPuppiWeightProducer::L1MultPuppiWeightProducer(const edm::ParameterSet& iConfig) :
  pfToken_(consumes<std::vector<l1t::PFCandidate>>(iConfig.getParameter<edm::InputTag>("pf"))),
  pupToken_(consumes<std::vector<l1t::PFCandidate>>(iConfig.getParameter<edm::InputTag>("pup"))),
  pupToken1_(consumes<std::vector<l1t::PFCandidate>>(iConfig.getParameter<edm::InputTag>("pup1"))),
  pupToken2_(consumes<std::vector<l1t::PFCandidate>>(iConfig.getParameter<edm::InputTag>("pup2"))),
  pupToken3_(consumes<std::vector<l1t::PFCandidate>>(iConfig.getParameter<edm::InputTag>("pup3"))),
  pupToken4_(consumes<std::vector<l1t::PFCandidate>>(iConfig.getParameter<edm::InputTag>("pup4"))),
  top_m_(iConfig.getParameter<unsigned int>("top_m"))
{
  produces<std::vector<l1t::PFCandidate>>("MultiPuppi");
}

L1MultPuppiWeightProducer::~L1MultPuppiWeightProducer() { }

int
L1MultPuppiWeightProducer::getMatch(int i, l1t::PFCandidate pfCand, std::vector<l1t::PFCandidate> pupCands) const
{
  float pfEta = pfCand.eta();
  float pfPhi = pfCand.phi();
  float pupEta = 999;
  float pupPhi = 999;

  bool found = false;


  for (size_t j = 0; j < pupCands.size(); ++j) {
    l1t::PFCandidate pupCand = (pupCands)[j];

    pupEta = pupCand.eta();
    pupPhi = pupCand.phi();

    bool kindCond = pupCand.id() == pfCand.id();

    found = (pfEta == pupEta and pfPhi == pupPhi and kindCond);

    if (found) {
      return j;
    }
  }

  return -1;
}

// ------------ method called for each event  ------------
void
L1MultPuppiWeightProducer::produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const
{
  edm::Handle<std::vector<l1t::PFCandidate>> pfParticles;
  iEvent.getByToken(pfToken_, pfParticles);

  edm::Handle<std::vector<l1t::PFCandidate>> pupParticles;
  iEvent.getByToken(pupToken_, pupParticles);

  edm::Handle<std::vector<l1t::PFCandidate>> pupParticles1;
  iEvent.getByToken(pupToken1_, pupParticles1);

  edm::Handle<std::vector<l1t::PFCandidate>> pupParticles2;
  iEvent.getByToken(pupToken2_, pupParticles2);

  edm::Handle<std::vector<l1t::PFCandidate>> pupParticles3;
  iEvent.getByToken(pupToken3_, pupParticles3);

  edm::Handle<std::vector<l1t::PFCandidate>> pupParticles4;
  iEvent.getByToken(pupToken4_, pupParticles4);

  std::vector<l1t::PFCandidate> cands;

  for (size_t i = 0; i < pfParticles->size(); ++i)
  {
    vector<l1t::PFCandidate> cand_list;
    vector<vector<l1t::PFCandidate>> pupCands;

    pupCands.emplace_back(*pupParticles);
    pupCands.emplace_back(*pupParticles1);
    pupCands.emplace_back(*pupParticles2);
    pupCands.emplace_back(*pupParticles3);
    pupCands.emplace_back(*pupParticles4);

    size_t max_j_num = top_m_;
    if (top_m_ > pupCands.size()) max_j_num = pupCands.size();

    int cand_i = -1;
    for (size_t j = 0; j < max_j_num; ++j) {
      cand_i = getMatch(i, (*pfParticles)[i], pupCands[j]);
      // cand_i will be -1 if no matching PUP candidate is found
      if (cand_i >= 0) cand_list.emplace_back(pupCands[j][cand_i]);
    }

    // store pf particle and weights
    double max_weight = 0;
    double max_j = -1;
    for (size_t j = 0; j < cand_list.size(); ++j) {
      if (cand_list[j].puppiWeight() > max_weight) {
        max_weight = cand_list[j].puppiWeight();
        max_j = j;
      }
    }

    if (max_weight > 0) {
      cands.emplace_back(cand_list[max_j]);
    }
  }

  iEvent.put(std::move(std::make_unique<std::vector<l1t::PFCandidate>>(cands)), "MultiPuppi");

  // should be able to check if all of the puppi candidates were matched
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1MultPuppiWeightProducer);
