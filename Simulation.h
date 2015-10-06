#ifndef _simulation_
#define _simulation_

#include "Track.h"
#include "Matrix.h"
#include "Propagation.h"
#include "Geometry.h"

void setupTrackByToyMC(SVector3& pos, SVector3& mom, SMatrixSym66& covtrk, HitVec& hits, MCHitInfoVec& initialhitinfo, unsigned int itrack, int& charge, const Geometry&, TSVec& initTSs);
#endif
