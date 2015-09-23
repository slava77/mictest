#ifndef _kalmanutils_mplex_
#define _kalmanutils_mplex_

#include "Track.h"
#include "Matrix.h"

void updateParametersMPlex(const MPlexLS &psErr,  const MPlexLV& psPar, const MPlexQI &inChg,
                           const MPlexHS &msErr,  const MPlexHV& msPar,
			   const MPlexQF &hitsRl, const MPlexQF& hitsXi,
			         MPlexLS &outErr,       MPlexLV& outPar);

void computeChi2MPlex(const MPlexLS &psErr,  const MPlexLV& psPar, const MPlexQI &inChg,
		      const MPlexHS &msErr,  const MPlexHV& msPar,
		      const MPlexQF &hitsRl, const MPlexQF& hitsXi,
                            MPlexQF& outChi2);

#endif
