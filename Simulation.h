#ifndef _simulation_
#define _simulation_

#include "Track.h"
#include "Matrix.h"
#include "Propagation.h"
#include "Geometry.h"

void setupTrackByToyMC(SVector3& pos, SVector3& mom, SMatrixSym66& covtrk, HitVec& hits, unsigned int itrack, 
		       int& charge, float pt, const Geometry&, HitVec& initHits, MCHitInfoVec& initialhitinfo);

void setupTrackFromTextFile(SVector3& pos, SVector3& mom, SMatrixSym66& covtrk, HitVec& hits, unsigned int itrack, 
			    int& charge, float pt, const Geometry&, HitVec& initHits, MCHitInfoVec& initialhitinfo);

#endif
