#ifndef _buildtest_mplex_
#define _buildtest_mplex_

#include "Event.h"
#include "Track.h"

#ifdef TEST_CLONE_ENGINE
#include "CandCloner.h"
#endif

/*
void   generateTracks(std::vector<Track>& simtracks, int Ntracks);
void   make_validation_tree(const char         *fname,
                            std::vector<Track> &simtracks,
                            std::vector<Track> &rectracks);
*/

double runBuildingTest(std::vector<Track>& simtracks/*, std::vector<Track>& rectracks*/);
double runBuildingTestBestHit(std::vector<Track>& simtracks/*, std::vector<Track>& rectracks*/);

void buildTestParallel(std::vector<Track>& evt_seeds,std::vector<Track>& evt_track_candidates,
		       std::vector<std::vector<Hit> >& evt_lay_hits,std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,
		       const int& nhits_per_seed,const unsigned int& maxCand,const float& chi2Cut,const float& nSigma,const float& minDPhi);

void processCandidates(std::pair<Track, TrackState>& cand,std::vector<std::pair<Track, TrackState> >& tmp_candidates,
		       unsigned int ilay,std::vector<std::vector<Hit> >& evt_lay_hits,std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,
		       const int& nhits_per_seed,const unsigned int& maxCand,const float& chi2Cut,const float& nSigma,const float& minDPhi);

#ifdef TEST_CLONE_ENGINE
double runBuildingTestPlex(std::vector<Track>& simtracks/*, std::vector<Track>& rectracks*/, CandCloner& cloner);
#else
double runBuildingTestPlex(std::vector<Track>& simtracks/*, std::vector<Track>& rectracks*/);
#endif
double runBuildingTestPlexOld(std::vector<Track>& simtracks/*, std::vector<Track>& rectracks*/);
double runBuildingTestPlexBestHit(std::vector<Track>& simtracks/*, std::vector<Track>& rectracks*/);

#endif
