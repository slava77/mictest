#ifndef _hit_
#define _hit_

#include <cmath>

#include "Matrix.h"
#include <atomic>

namespace Config
{
  static const float    PI = 3.14159265358979323846;
  static const float TwoPI = 6.28318530717958647692;
  
  static constexpr const int   nPhiPart   = 1260; // 63;
  static constexpr const float nPhiFactor = nPhiPart / TwoPI;
  static constexpr const int   nEtaPart   = 11; // 10;

  static const int   nEtaBin   = 2*nEtaPart - 1;
  static const float fEtaDet   = 1;
  static const float fEtaFull  = 2 * fEtaDet;
  static const float lEtaPart  = fEtaFull/float(nEtaPart);
  static const float lEtaBin   = lEtaPart/2.;

  static const float fEtaOffB1 = fEtaDet;
  static const float fEtaFacB1 = nEtaPart / fEtaFull;
  static const float fEtaOffB2 = fEtaDet - fEtaFull / (2 * nEtaPart);
  static const float fEtaFacB2 = (nEtaPart - 1) / (fEtaFull - fEtaFull / nEtaPart);

  // This is for extra bins narrower ... thinking about this some more it
  // seems it would be even better to have two more exta bins, hanging off at
  // both ends.
  //
  // Anyway, it doesn't matter ... as with wide vertex region this eta binning
  // won't make much sense -- will have to be done differently for different
  // track orgin hypotheses. In about a year or so.

  inline int getEtaBin(float eta)
  {

    //in this case we are out of bounds
    if (fabs(eta)>fEtaDet) return -1;

    //first and last bin have extra width
    if (eta<(lEtaBin-fEtaDet)) return 0;
    if (eta>(fEtaDet-lEtaBin)) return nEtaBin-1;

    //now we can treat all bins as if they had same size
    return int( (eta+fEtaDet-lEtaBin/2.)/lEtaBin );

  }

  inline int getBothEtaBins(float eta, int& b1, int& b2)
  {
    b1 = b2 = -1;

    if (eta < -fEtaDet || eta > fEtaDet)
    {
      return 0;
    }

    int b1p = std::floor((eta + fEtaOffB1) * fEtaFacB1);
    int b2p = std::floor((eta + fEtaOffB2) * fEtaFacB2);

    // printf("b1' = %d   b2' = %d\n", b1p, b2p);

    int cnt = 0;
    if (b1p >= 0 && b1p < nEtaPart)
    {
      b1 = 2 * b1p;
      ++cnt;
    }
    if (b2p >= 0 && b2p < nEtaPart - 1)
    {
      b2 = 2 * b2p + 1;
      ++cnt;
    }

    // printf("b1  = %d   b2  = %d\n", b1, b2);

    return cnt;
  }
};

class HitsOnLayer
{
  
public:
};


inline int getPhiPartition(float phi)
{
  //assume phi is between -PI and PI
  //  if (!(fabs(phi)<TMath::Pi())) std::cout << "anomalous phi=" << phi << std::endl;
  float phiPlusPi  = phi + Config::PI;
  int bin = phiPlusPi * Config::nPhiFactor;
  return bin;
}

inline int getEtaPartition(float eta, float etaDet)
{
  float etaPlusEtaDet  = eta + etaDet;
  float twiceEtaDet    = 2.0 * etaDet;
  int bin     = (etaPlusEtaDet * Config::nEtaPart) / twiceEtaDet;  // ten bins for now ... update if changed in Event.cc
  return bin;
}

inline float getPhi(float x, float y)
{
  return std::atan2(y,x); 
}

inline float getEta(float x, float y, float z)
{
  float theta = atan2( std::sqrt(x*x+y*y), z );
  return -1. * log( tan(theta/2.) );
}

inline float getHypot(float x, float y)
{
  return sqrtf(x*x + y*y);
}


struct MCHitInfo
{
  MCHitInfo() : mcHitID_(mcHitIDCounter_++) {}
  MCHitInfo(int track, int layer, int ithlayerhit)
    : mcTrackID_(track), layer_(layer), ithLayerHit_(ithlayerhit), mcHitID_(++mcHitIDCounter_) {}

  int mcTrackID_;
  int layer_;
  int ithLayerHit_;
  int mcHitID_;

  static std::atomic<unsigned int> mcHitIDCounter_;
};

struct MeasurementState
{
  SVector3 parameters() const { return SVector3(pos[0],pos[1],pos[2]); }
  SMatrixSym33 errors() const { 
    SMatrixSym33 result;
    for (int i=0;i<6;++i) result.Array()[i]=err[i];
    return result; 
  }
public:
  float pos[3];
  float err[6];
};

class Hit
{
public:
  Hit(){}

  Hit(const MeasurementState& state) : state_(state) {}

  Hit(const SVector3& position, const SMatrixSym33& error, int mcHitTrkID)
  {
    state_.pos[0]=position[0];
    state_.pos[1]=position[1];
    state_.pos[2]=position[2];
    state_.err[0]=error.Array()[0];
    state_.err[1]=error.Array()[1];
    state_.err[2]=error.Array()[2];
    state_.err[3]=error.Array()[3];
    state_.err[4]=error.Array()[4];
    state_.err[5]=error.Array()[5];
    mcHitTrkID_ = mcHitTrkID;
  }

  ~Hit(){}

  const SVector3     position()   const {return state_.parameters();}
  const SVector3     parameters() const {return state_.parameters();}
  const SMatrixSym33 error()      const {return state_.errors();}

  // Non-const versions needed for CopyOut of Matriplex.
  SVector3     parameters_nc() {return state_.parameters();}
  SMatrixSym33 error_nc()      {return state_.errors();}

  float r() const {
    return sqrt(state_.pos[0]*state_.pos[0] +
                state_.pos[1]*state_.pos[1]);
  }
  float x() const {
    return state_.pos[0];
  }
  float y() const {
    return state_.pos[1];
  }
  float z() const {
    return state_.pos[2];
  }
  float phi() const {
    return getPhi(state_.pos[0], state_.pos[1]);
  }
  float eta() const {
    return getEta(state_.pos[0], state_.pos[1], state_.pos[2]);
  }

  int mcHitTrkID() const { return mcHitTrkID_; }
  int mcTrackID() const { return mcHitTrkID_ / 10; }
  int mcHitID() const { return mcHitTrkID_ % 10; }

  const MeasurementState& measurementState() const {
    return state_;
  }

private:
  MeasurementState state_;
  // unique hit index per event, defined as (10*simtrack_idex + 1*hit_idx). 
  // the factor 10* is ok since we have 10 hits per track, it can become 100* if we allow for more. 
  // we will need a different indexing when we'll consider htis from real detector
  int mcHitTrkID_;
};

typedef std::vector<Hit> HitVec;
typedef std::vector<MCHitInfo> MCHitInfoVec;

#endif
