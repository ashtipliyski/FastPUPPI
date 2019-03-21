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

  float getWeight(int i, l1t::PFCandidate pfCand, std::vector<l1t::PFCandidate> pupCands) const;
  l1t::PFCandidate collapseWeights(l1t::PFCandidate pfCand, std::vector<double> weights);

  // edm::EDGetTokenT<reco::CandidateView> gen_;
  // StringCutObjectSelector<reco::Candidate> sel_;

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

float
L1MultPuppiWeightProducer::getWeight(int i, l1t::PFCandidate pfCand, std::vector<l1t::PFCandidate> pupCands) const
{
  float weight = 0;

  float pfEta = pfCand.eta();
  float pfPhi = pfCand.phi();
  float pupEta = 999;
  float pupPhi = 999;

  bool found = false;


  for (size_t j = 0; j < pupCands.size(); ++j) {
    pupEta = (pupCands)[j].eta();
    pupPhi = (pupCands)[j].phi();

    found = (pfEta == pupEta and pfPhi == pupPhi);

    if (found) {
      weight = pupCands[j].puppiWeight();
      break;
    }
  }

  if (!found) {
    weight = 0;
  }

  return weight;
}

l1t::PFCandidate collapseWeights(l1t::PFCandidate pfCand, std::vector<double> weights)
{
  // current implementation picks the max weight and attaches it to candidate
  double weight = 0;
  for (double w : weights) weight = std::max(weight, w);

  l1t::PFCandidate cand(pfCand);
  cand.setPuppiWeight(weight);

  return cand;
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
    vector<double> weights;

    weights.emplace_back(getWeight(i, (*pfParticles)[i], *pupParticles));
    weights.emplace_back(getWeight(i, (*pfParticles)[i], *pupParticles1));
    weights.emplace_back(getWeight(i, (*pfParticles)[i], *pupParticles2));
    weights.emplace_back(getWeight(i, (*pfParticles)[i], *pupParticles3));
    weights.emplace_back(getWeight(i, (*pfParticles)[i], *pupParticles4));

    double weight = 0;
    for (const auto & w : weights) weight = std::max(weight, w);

    l1t::PFCandidate cand((*pfParticles)[i]);
    cand.setPuppiWeight(weight);
    cands.emplace_back(cand);

    // store pf particle and weights
  }

  iEvent.put(std::move(std::make_unique<std::vector<l1t::PFCandidate>>(cands)), "MultiPuppi");

  // should be able to check if all of the puppi candidates were matched
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1MultPuppiWeightProducer);
