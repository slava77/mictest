#include "buildtest.h"

#include "KalmanUtils.h"
#include "Propagation.h"
#include "Simulation.h"
#include "Event.h"
#include "Debug.h"

#include <cmath>
#include <iostream>

//typedef std::pair<Track, TrackState> cand_t;
typedef Track cand_t;
typedef TrackVec::const_iterator TrkIter;

#ifndef TBB
typedef std::vector<cand_t> candvec;
#else
#include "tbb/tbb.h"
#include <mutex>
// concurrent_vector is only needed if we parallelize the candidate loops;
// not needed if we only parallelize over seeds.
//typedef tbb::concurrent_vector<cand_t> candvec;
typedef std::vector<cand_t> candvec;
static std::mutex evtlock;
#endif
typedef candvec::const_iterator canditer;

void extendCandidate(const Event& ev, const cand_t& cand, candvec& tmp_candidates, unsigned int ilay, bool debug);

#ifdef ETASEG
inline float normalizedEta(float eta) {
  static float const ETA_DET = 2.0;

  if (eta < -ETA_DET ) eta = -ETA_DET+.00001;
  if (eta >  ETA_DET ) eta =  ETA_DET-.00001;
  return eta;
}
#endif

static bool sortByHitsChi2(const cand_t& cand1, const cand_t& cand2)
{
  if (cand1.nFoundHits()==cand2.nFoundHits()) return cand1.chi2()<cand2.chi2();
  return cand1.nFoundHits()>cand2.nFoundHits();
}

void processCandidates(Event& ev, const Track& seed, candvec& candidates, unsigned int ilay, const bool debug)
{
  auto& evt_track_candidates(ev.candidateTracks_);

  dprint("processing seed # " << seed.SimTrackIDInfo().first << " par=" << seed.parameters() << " candidates=" << candidates.size());

  candvec tmp_candidates;
  tmp_candidates.reserve(3*candidates.size()/2);

  if (candidates.size() > 0) {
    //loop over running candidates
    for (auto&& cand : candidates) {
      extendCandidate(ev, cand, tmp_candidates, ilay, debug);
    }

    ev.validation_.fillBuildHists(ilay, tmp_candidates.size(), candidates.size());

    if (tmp_candidates.size()>Config::maxCand) {
      dprint("huge size=" << tmp_candidates.size() << " keeping best "<< Config::maxCand << " only");
      std::partial_sort(tmp_candidates.begin(),tmp_candidates.begin()+(Config::maxCand-1),tmp_candidates.end(),sortByHitsChi2);
      tmp_candidates.resize(Config::maxCand); // thread local, so ok not thread safe
    } else if (tmp_candidates.size()==0) {
      dprint("no more candidates, saving best");
      // save the best candidate from the previous iteration and then swap in
      // the empty new candidate list; seed will be skipped on future iterations
      auto&& best = std::min_element(candidates.begin(),candidates.end(),sortByHitsChi2);
#ifdef TBB
      std::lock_guard<std::mutex> evtguard(evtlock); // should be rare
#endif
      evt_track_candidates.push_back(*best);
    }
    dprint("swapping with size=" << tmp_candidates.size());
    candidates.swap(tmp_candidates);
    tmp_candidates.clear();
  }
}

void buildTracksBySeeds(Event& ev)
{
  auto& evt_track_candidates(ev.candidateTracks_);
  const auto& evt_lay_hits(ev.layerHits_);
  const auto& evt_seeds(ev.seedTracks_);
  bool debug(false);

  std::vector<candvec> track_candidates(evt_seeds.size());
  for (auto iseed = 0U; iseed < evt_seeds.size(); iseed++) {
    const auto& seed(evt_seeds[iseed]);
    track_candidates[iseed].push_back(seed);
  }

#ifdef TBB
  //loop over seeds
  parallel_for( tbb::blocked_range<size_t>(0, evt_seeds.size()), 
      [&](const tbb::blocked_range<size_t>& seediter) {
    for (auto iseed = seediter.begin(); iseed != seediter.end(); ++iseed) {
      const auto& seed(evt_seeds[iseed]);
#else
    for (auto iseed = 0U; iseed != evt_seeds.size(); ++iseed) {
      const auto& seed(evt_seeds[iseed]);
#endif
      dprint("processing seed # " << seed.SimTrackIDInfo().first << " par=" << seed.parameters());
      TrackState seed_state = seed.state();
      //seed_state.errors *= 0.01;//otherwise combinatorics explode!!!
      //should consider more than 1 candidate...
      auto&& candidates(track_candidates[iseed]);
      for (unsigned int ilay=Config::nlayers_per_seed;ilay<evt_lay_hits.size();++ilay) {//loop over layers, starting from after the seed
        dprint("going to layer #" << ilay << " with N cands=" << track_candidates.size());
        processCandidates(ev, seed, candidates, ilay, debug);
      }
      //end of layer loop
    }//end of process seeds loop
#ifdef TBB
  });
#endif
  for (const auto& cand : track_candidates) {
    if (cand.size()>0) {
      // only save one track candidate per seed, one with lowest chi2
      //std::partial_sort(cand.begin(),cand.begin()+1,cand.end(),sortByHitsChi2);
      auto&& best = std::min_element(cand.begin(),cand.end(),sortByHitsChi2);
      evt_track_candidates.push_back(*best);
    }
  }
}

void buildTracksByLayers(Event& ev)
{
  auto& evt_track_candidates(ev.candidateTracks_);
  const auto& evt_lay_hits(ev.layerHits_);
  const auto& evt_seeds(ev.seedTracks_);
  bool debug(false);

  std::vector<candvec> track_candidates(evt_seeds.size());
  for (auto iseed = 0U; iseed < evt_seeds.size(); iseed++) {
    const auto& seed(evt_seeds[iseed]);
    track_candidates[iseed].push_back(seed);
  }

  //loop over layers, starting from after the seed
  for (auto ilay = Config::nlayers_per_seed; ilay < evt_lay_hits.size(); ++ilay) {
    dprint("going to layer #" << ilay << " with N cands = " << track_candidates.size());

#ifdef TBB
    //loop over seeds
    parallel_for( tbb::blocked_range<size_t>(0, evt_seeds.size()), 
        [&](const tbb::blocked_range<size_t>& seediter) {
      for (auto iseed = seediter.begin(); iseed != seediter.end(); ++iseed) {
        const auto& seed(evt_seeds[iseed]);
        auto&& candidates(track_candidates[iseed]);
        processCandidates(ev, seed, candidates, ilay, debug);
      }
    }); //end of process seeds loop
#else
    //process seeds
    for (auto iseed = 0U; iseed != evt_seeds.size(); ++iseed) {
      const auto& seed(evt_seeds[iseed]);
      auto&& candidates(track_candidates[iseed]);
      processCandidates(ev, seed, candidates, ilay, debug);
    }
#endif
  } //end of layer loop

  //std::lock_guard<std::mutex> evtguard(evtlock);
  for (const auto& cand : track_candidates) {
    if (cand.size()>0) {
      // only save one track candidate per seed, one with lowest chi2
      //std::partial_sort(cand.begin(),cand.begin()+1,cand.end(),sortByHitsChi2);
      auto&& best = std::min_element(cand.begin(),cand.end(),sortByHitsChi2);
      evt_track_candidates.push_back(*best);
    }
  }
}

void extendCandidate(const Event& ev, const cand_t& cand, candvec& tmp_candidates, unsigned int ilayer,
                     bool debug)
{
  const Track& tkcand = cand;
  const TrackState& updatedState = cand.state();
  const auto& evt_lay_hits(ev.layerHits_);
  const auto& segmentMap(ev.segmentMap_);
  //  debug = true;

  dprint("processing candidate with nHits=" << tkcand.nHits());
#ifdef LINEARINTERP
  TrackState propState = propagateHelixToR(updatedState,ev.geom_.Radius(ilayer));
#else
#ifdef TBB
#error "Invalid combination of options (thread safety)"
#endif
  TrackState propState = propagateHelixToLayer(updatedState,ilayer,ev.geom_);
#endif // LINEARINTERP
#ifdef CHECKSTATEVALID
  if (!propState.valid) {
    return;
  }
#endif
  const float predx = propState.parameters.At(0);
  const float predy = propState.parameters.At(1);
  const float predz = propState.parameters.At(2);
  const float px2py2 = predx*predx+predy*predy; // predicted radius^2
#ifdef DEBUG
  if (debug) {
    std::cout << "propState at hit#" << ilayer << " r/phi/z : " << sqrt(pow(predx,2)+pow(predy,2)) << " "
              << std::atan2(predy,predx) << " " << predz << std::endl;
    dumpMatrix(propState.errors);
  }
#endif

#ifdef ETASEG  
  const float eta = getEta(predx,predy,predz);
  const float pz2 = predz*predz;
  const float detadx = -predx/(px2py2*sqrt(1+(px2py2/pz2)));
  const float detady = -predy/(px2py2*sqrt(1+(px2py2/pz2)));
  const float detadz = 1.0/(predz*sqrt(1+(px2py2/pz2)));
  const float deta2  = detadx*detadx*(propState.errors.At(0,0)) +
    detady*detady*(propState.errors.At(1,1)) +
    detadz*detadz*(propState.errors.At(2,2)) +
    2*detadx*detady*(propState.errors.At(0,1)) +
    2*detadx*detadz*(propState.errors.At(0,2)) +
    2*detady*detadz*(propState.errors.At(1,2));
  const float deta   = sqrt(std::abs(deta2));  
  const float nSigmaDeta = std::min(std::max(Config::nSigma*deta,(float)0.0),float( 1.0)); // something to tune -- minDEta = 0.0

  //for now as well --> eta boundary!!!
  const float detaMinus  = normalizedEta(eta-nSigmaDeta);
  const float detaPlus   = normalizedEta(eta+nSigmaDeta);
  
  // for now
  const float etaDet = 2.0;

  const auto etaBinMinus = getEtaPartition(detaMinus,etaDet);
  const auto etaBinPlus  = getEtaPartition(detaPlus,etaDet);

  dprint("eta: " << eta << " etaBinMinus: " << etaBinMinus << " etaBinPlus: " << etaBinPlus << " deta2: " << deta2);
#else // just assign the etaBin boundaries to keep code without 10k ifdefs
  const auto etaBinMinus = 0U;
  const auto etaBinPlus  = 0U;
#endif

  const float phi = getPhi(predx,predy); //std::atan2(predy,predx); 
  const float dphidx = -predy/px2py2;
  const float dphidy =  predx/px2py2;
  const float dphi2  = dphidx*dphidx*(propState.errors.At(0,0)) +
    dphidy*dphidy*(propState.errors.At(1,1)) +
    2*dphidx*dphidy*(propState.errors.At(0,1));
  
  const float dphi   =  sqrt(std::abs(dphi2));//how come I get negative squared errors sometimes?
  const float nSigmaDphi = std::min(std::max(Config::nSigma*dphi,(float) Config::minDPhi), (float) M_PI);
  
  const float dphiMinus = normalizedPhi(phi-nSigmaDphi);
  const float dphiPlus  = normalizedPhi(phi+nSigmaDphi);
  
  const auto phiBinMinus = getPhiPartition(dphiMinus);
  const auto phiBinPlus  = getPhiPartition(dphiPlus);

  dprint("phi: " << phi << " phiBinMinus: " <<  phiBinMinus << " phiBinPlus: " << phiBinPlus << "  dphi2: " << dphi2);
  
  for (unsigned int ieta = etaBinMinus; ieta <= etaBinPlus; ++ieta){
    
    const BinInfo binInfoMinus = segmentMap[ilayer][ieta][int(phiBinMinus)];
    const BinInfo binInfoPlus  = segmentMap[ilayer][ieta][int(phiBinPlus)];
    
    std::vector<unsigned int> cand_hit_idx;
    std::vector<unsigned int>::iterator index_iter; // iterator for vector
    
    // Branch here from wrapping
    if (phiBinMinus<=phiBinPlus){
      const auto firstIndex = binInfoMinus.first;
      const auto maxIndex   = binInfoPlus.first+binInfoPlus.second;

      for (auto ihit  = firstIndex; ihit < maxIndex; ++ihit){
        cand_hit_idx.push_back(ihit);
      }
    } 
    else { // loop wrap around end of array for phiBinMinus > phiBinPlus, for dPhiMinus < 0 or dPhiPlus > 0 at initialization
      const auto firstIndex = binInfoMinus.first;
      const auto etaBinSize = segmentMap[ilayer][ieta][Config::nPhiPart-1].first+segmentMap[ilayer][ieta][Config::nPhiPart-1].second;

      for (auto ihit  = firstIndex; ihit < etaBinSize; ++ihit){
        cand_hit_idx.push_back(ihit);
      }

      const auto etaBinStart= segmentMap[ilayer][ieta][0].first;
      const auto maxIndex   = binInfoPlus.first+binInfoPlus.second;

      for (unsigned int ihit  = etaBinStart; ihit < maxIndex; ++ihit){
        cand_hit_idx.push_back(ihit);
      }
    }
  
#ifdef LINEARINTERP
    const float minR = ev.geom_.Radius(ilayer);
    float maxR = minR;
    for(index_iter = cand_hit_idx.begin(); index_iter != cand_hit_idx.end(); ++index_iter){
      const float candR = evt_lay_hits[ilayer][*index_iter].r();
      if (candR > maxR) maxR = candR;
    }
    const float deltaR = maxR - minR;
    dprint("min, max: " << minR << ", " << maxR);
    const TrackState propStateMin = propState;
    const TrackState propStateMax = propagateHelixToR(updatedState,maxR);
#ifdef CHECKSTATEVALID
    if (!propStateMax.valid) {
      return;
    }
#endif
#endif
    
    for(index_iter = cand_hit_idx.begin(); index_iter != cand_hit_idx.end(); ++index_iter){
      const Hit hitCand = evt_lay_hits[ilayer][*index_iter];
      const MeasurementState hitMeas = hitCand.measurementState();

#ifdef LINEARINTERP
      const float ratio = (hitCand.r() - minR)/deltaR;
      propState.parameters = (1.0-ratio)*propStateMin.parameters + ratio*propStateMax.parameters;
      dprint(std::endl << ratio << std::endl << propStateMin.parameters << std::endl << propState.parameters << std::endl
                       << propStateMax.parameters << std::endl << propStateMax.parameters - propStateMin.parameters
                       << std::endl << std::endl << hitMeas.parameters);
#endif
      const float chi2 = computeChi2(propState,hitMeas);
    
      if ((chi2<Config::chi2Cut)&&(chi2>0.)) {//fixme 
        dprint("found hit with index: " << *index_iter << " chi2=" << chi2);
        const TrackState tmpUpdatedState = updateParameters(propState, hitMeas);
        Track tmpCand = tkcand.clone();//(tmpUpdatedState,tkcand.hitsVector(),tkcand.chi2()); //= tkcand.clone();
        tmpCand.addHitIdx(*index_iter,chi2);
        tmp_candidates.push_back(tmpCand);
      }
    }//end of consider hits on layer loop
    
  }//end of eta loop

  //add also the candidate for no hit found
  if (tkcand.nFoundHits()==ilayer) {//only if this is the first missing hit
    dprint("adding candidate with no hit");
    tmp_candidates.push_back(tkcand);
  }
}
