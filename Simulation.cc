#include <cmath>

#include "Simulation.h"
#include "Debug.h"

//#define SOLID_SMEAR
#define SCATTER_XYZ

void setupTrackByToyMC(SVector3& pos, SVector3& mom, SMatrixSym66& covtrk, HitVec& hits, unsigned int itrack,
                       int& charge, float pt, const Geometry& geom, HitVec& initHits, MCHitInfoVec& initialhitinfo)
{
#ifdef DEBUG
  bool debug = false;
#endif
  //assume beam spot width 1mm in xy and 1cm in z
  pos=SVector3(0.1*g_gaus(g_gen), 0.1*g_gaus(g_gen), 1.0*g_gaus(g_gen));

  dprint("pos x=" << pos[0] << " y=" << pos[1] << " z=" << pos[2]);

  if (charge==0) {
    if (g_unif(g_gen) > 0.5) charge = -1;
    else charge = 1;
  }

  //float phi = 0.5*TMath::Pi()*(1-g_unif(g_gen)); // make an angle between 0 and pi/2 //fixme
  float phi = 0.5*TMath::Pi()*g_unif(g_gen); // make an angle between 0 and pi/2
  float px = pt * cos(phi);
  float py = pt * sin(phi);

  if (g_unif(g_gen)>0.5) px*=-1.;
  if (g_unif(g_gen)>0.5) py*=-1.;
  float pz = pt*(2.3*(g_unif(g_gen)-0.5));//so that we have -1<eta<1

  mom=SVector3(px,py,pz);
  covtrk=ROOT::Math::SMatrixIdentity();
  //initial covariance can be tricky
  for (unsigned int r=0; r<6; ++r) {
    for (unsigned int c=0; c<6; ++c) {
      if (r==c) {
	if (r<3) covtrk(r,c)=pow(1.0*pos[r],2);//100% uncertainty on position
	else covtrk(r,c)=pow(1.0*mom[r-3],2);  //100% uncertainty on momentum
      } else {
        covtrk(r,c)=0.;                   //no covariance
      }
    }
  }

  dprint("track with px: " << px << " py: " << py << " pz: " << pz << " pt: " << sqrt(px*px+py*py) << " p: " << sqrt(px*px+py*py+pz*pz) << std::endl);

  const float hitposerrXY = 0.01;//assume 100mum uncertainty in xy coordinate
  const float hitposerrZ = 0.1;//assume 1mm uncertainty in z coordinate
  const float hitposerrR = hitposerrXY/10.;

  const float varXY  = hitposerrXY*hitposerrXY;
  const float varZ   = hitposerrZ*hitposerrZ;

  TrackState initState;
  initState.parameters=SVector6(pos[0],pos[1],pos[2],mom[0],mom[1],mom[2]);
  initState.errors=covtrk;
  initState.charge=charge;

  TrackState tmpState = initState;

  // useful info for loopers/overlaps

  unsigned int nLayers = geom.CountLayers();
  unsigned int layer_counts[nLayers];
  for (unsigned int ilayer=0;ilayer<nLayers;++ilayer){
    layer_counts[ilayer]=0;
  }

  unsigned int nTotHit = nLayers; // can tune this number!
  // to include loopers, and would rather add a break on the code if layer ten exceeded
  // if block BREAK if hit.Layer == theGeom->CountLayers() 
  // else --> if (NMAX TO LOOPER (maybe same as 10?) {break;} else {continue;}
  
  unsigned int simLayer = 0;

  hits.reserve(nTotHit);
  initHits.reserve(nTotHit);

  for (unsigned int ihit=0;ihit<nTotHit;++ihit) {  // go to first layer in radius using propagation.h
    //TrackState propState = propagateHelixToR(tmpState,4.*float(ihit+1));//radius of 4*ihit
    auto propState = propagateHelixToNextSolid(tmpState,geom);

    float initX   = propState.parameters.At(0);
    float initY   = propState.parameters.At(1);
    float initZ   = propState.parameters.At(2);
    float initPhi = atan2(initY,initX);
    float initRad = sqrt(initX*initX+initY*initY);

    UVector3 init_point(initX,initY,initZ);
    simLayer = geom.LayerIndex(init_point);

#ifdef SCATTERING
    // PW START
    // multiple scattering. Follow discussion in PDG Ch 32.3
    // this generates a scattering distance yplane and an angle theta_plane in a coordinate
    // system normal to the initial direction of the incident particle ...
    // question: do I need to generate two scattering angles in two planes (theta_plane) and add those or 
    // can I do what I did, generate one (theta_space) and rotate it by another random numbers
    const float z1 = g_gaus(g_gen);
    const float z2 = g_gaus(g_gen);
    const float phismear = g_unif(g_gen)*TMath::TwoPi(); // random rotation of scattering plane
    const float X0 = 9.370; // cm, from http://pdg.lbl.gov/2014/AtomicNuclearProperties/HTML/silicon_Si.html
    //const float X0 = 0.5612; // cm, for Pb
    float x = 0.1; // cm  -assumes radial impact. This is bigger than what we have in main
    // will update for tilt down a few lines
    const float p = sqrt(propState.parameters.At(3)*propState.parameters.At(3)+
                         propState.parameters.At(4)*propState.parameters.At(4)+
                         propState.parameters.At(5)*propState.parameters.At(5));
    UVector3 pvec(propState.parameters.At(3)/p, propState.parameters.At(4)/p, propState.parameters.At(5)/p);

    // now we need to transform to the right coordinate system
    // in this coordinate system z -> z'=z i.e., unchanged (this is a freedom I have since the 
    // definition of the plane defined by nhat is invariant under rotations about nhat
    // x, y change
    // so the transformation is a rotation about z and a translation
    // if the normal vector is given by n == x' = (x1,y1,0) [NB no component along z here]
    // then the 
    // y axis in the new coordinate system is given by y' = z' x x' = (-y1,x1,0)
    const auto theInitSolid = geom.InsideWhat(init_point);
    if ( ! theInitSolid ) {
      std::cerr << __FILE__ << ":" << __LINE__ << ": failed to find solid." <<std::endl;
      std::cerr << "itrack = " << itrack << ", ihit = " << ihit << ", r = " << initRad << ", r*4cm = " << 4*ihit << ", phi = " << initPhi << std::endl;
      std::cerr << "initX = " << initX << ", initY = " << initY << ", initZ = " << initZ << std::endl;
      float pt = sqrt(propState.parameters[3]*propState.parameters[3]+
                      propState.parameters[4]*propState.parameters[4]);
      std::cerr << "pt = " << pt << ", pz = " << propState.parameters[5] << std::endl;

      continue;
    }
    UVector3 init_xprime; // normal on surface
    bool init_good = theInitSolid->Normal(init_point, init_xprime);
      
    if ( ! init_good ) {
      std::cerr << __FILE__ << ":" << __LINE__ << ": failed to find normal vector at " << init_point <<std::endl;
      break;
    }
    assert(std::abs(init_xprime[2])<1e-10); // in this geometry normal vector is in xy plane

    x /= std::abs(init_xprime.Dot(pvec)); // take care of dip angle
    const float betac = sqrt(p*p+(.135*.135))/(p*p); // m=130 MeV, pion
    const float theta_0 = 0.0136/(betac*p)*sqrt(x/X0)*(1+0.038*log(x/X0));// eq 32.14
    const float y_plane = z1*x*theta_0/sqrt(12.)+ z2*x*theta_0/2.;
    const float theta_plane = z2*theta_0;
    const float theta_space = sqrt(2)*theta_plane;
    dprint("yplane, theta_space = " << y_plane << ", " << theta_space);

    UVector3 yprime(-init_xprime[1],init_xprime[0],0); // result of cross product with zhat
    //const double phi0 = atan2(xpime[1], init_xprime[0]);

#ifdef SCATTER_XYZ
    const float scatteredX = initX + y_plane *(-init_xprime[1]*cos(phismear)); 
    const float scatteredY = initY + y_plane *(+init_xprime[0]*cos(phismear));
    const float scatteredZ = initZ + y_plane *(           sin(phismear));
    const float scatteredPhi = atan2(scatteredY,scatteredX);
    const float scatteredRad = sqrt(scatteredX*scatteredX+scatteredY*scatteredY);
#else 
    const float scatteredX = initX;
    const float scatteredY = initY;
    const float scatteredZ = initZ;
    const float scatteredPhi = atan2(initY,initX);
    const float scatteredRad = sqrt(initX*initX+initY*initY);
#endif  // SCATTER_XYZ
    UVector3 pvecprime;

    const float v0 = sqrt(2+pow((pvec[1]+pvec[2])/pvec[0],2));
    const float v1 = sqrt(2-pow((pvec[1]-pvec[2]),2));
    const float a = pvec[0]; 
    const float b = pvec[1]; 
    const float c = pvec[2];

    auto sign = [] (const float a) { 
      if ( a > 0 ) 
        return 1;
      else if ( a < 0 ) 
        return -1;
      else
        return 0;
    };
                    
    pvecprime[0] = a*cos(theta_space) + ((b + c)*cos(phismear)*sin(theta_space))/(a*v0) + 
      (a*(b - c)*sin(theta_space)*sin(phismear))/(v1*sign(1 + (b - c)*c));
    pvecprime[1] = b*cos(theta_space) - (cos(phismear)*sin(theta_space))/v0 + 
      ((-1 + pow(b,2) - b*c)*sin(theta_space)*sin(phismear))/(v1*sign(1 + (b - c)*c));
    pvecprime[2] = c*cos(theta_space) - (cos(phismear)*sin(theta_space))/v0 + 
      (std::abs(1 + (b - c)*c)*sin(theta_space)*sin(phismear))/v1; 
    assert(pvecprime.Mag()<=1.0001);

    dprint("theta_space, phismear = " << theta_space << ", " << phismear << std::endl
        << "init_xprime = " << init_xprime << std::endl
        << "yprime = " << yprime << std::endl
        << "phat      = " << pvec << "\t" << pvec.Mag() << std::endl
        << "pvecprime = " << pvecprime << "\t" << pvecprime.Mag() << std::endl
        << "angle     = " << pvecprime.Dot(pvec) << "(" << cos(theta_space) << ")" << std::endl
        << "pt, before and after: " << pvec.Perp()*p << ", "<< pvecprime.Perp()*p);

    pvecprime.Normalize();

    // now update propstate with the new scattered results

    propState.parameters[0] = scatteredX;
    propState.parameters[1] = scatteredY;
    propState.parameters[2] = scatteredZ;

    propState.parameters[3] = pvecprime[0]*p;
    propState.parameters[4] = pvecprime[1]*p;
    propState.parameters[5] = pvecprime[2]*p;
    
    // PW END
    float hitZ    = hitposerrZ*g_gaus(g_gen)+scatteredZ;
    float hitPhi  = ((hitposerrXY/scatteredRad)*g_gaus(g_gen))+scatteredPhi;
    float hitRad  = (hitposerrR)*g_gaus(g_gen)+scatteredRad;

#else // SCATTERING
    float hitZ    = hitposerrZ*g_gaus(g_gen)+initZ;
    float hitPhi  = ((hitposerrXY/initRad)*g_gaus(g_gen))+initPhi;
    float hitRad  = (hitposerrR)*g_gaus(g_gen)+initRad;
#endif // SCATTERING

#ifdef SOLID_SMEAR
    UVector3 scattered_point(scatteredX,scatteredY,scatteredZ);
    const auto theScatteredSolid = geom.InsideWhat(scattered_point);
    if ( ! theScatteredSolid ) {
      std::cerr << __FILE__ << ":" << __LINE__ << ": failed to find solid AFTER scatter." << std::endl;
      std::cerr << "itrack = " << itrack << ", ihit = " << ihit << ", r = " << sqrt(scatteredX*scatteredX + scatteredY*scatteredY) << ", r*4cm = " << 4*ihit << ", phi = " << atan2(scatteredY,scatteredX) << std::endl;
      std::cerr << "initX = " << initX << ", initY = " << initY << ", initZ = " << initZ << std::endl;
      std::cerr << "scatteredX = " << scatteredX << ", scatteredY = " << scatteredY << ", scatteredZ = " << scatteredZ << std::endl << std::endl;       
      //    continue;
    }

    UVector3 scattered_xprime; // normal on surface
    bool scattered_good = theScatteredSolid->Normal(scattered_point, scattered_xprime);
      
    if ( ! scattered_good ) {
      std::cerr << __FILE__ << ":" << __LINE__ << ": failed to find normal vector at " << scattered_point <<std::endl;
    }
    assert(std::abs(scattered_xprime[2])<1e-10); // in this geometry normal vector is in xy plane
    
    // smear along scattered yprime (where yprime = z_prime x x_prime, zprime = (0,0,1)

    float xyres_smear = g_gaus(g_gen)*hitposerrXY;

    float hitX = scatteredX - (scattered_xprime[1] * xyres_smear); 
    float hitY = scatteredY + (scattered_xprime[0] * xyres_smear); 
#else
    float hitRad2 = hitRad*hitRad;
    float hitX    = hitRad*cos(hitPhi);
    float hitY    = hitRad*sin(hitPhi);

    float varPhi = varXY/hitRad2;
    float varR   = hitposerrR*hitposerrR;
#endif

#ifdef DEBUG
    if (debug) {
      UVector3 hit_point(hitX,hitY,hitZ);
      const auto theHitSolid = geom.InsideWhat(hit_point);
      if ( ! theHitSolid ) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": failed to find solid AFTER scatter+smear." << std::endl;
        std::cerr << "itrack = " << itrack << ", ihit = " << ihit << ", r = " << sqrt(hitX*hitX + hitY*hitY) << ", r*4cm = " << 4*ihit << ", phi = " << atan2(hitY,hitX) << std::endl;
        std::cerr << "initX = " << initX << ", initY = " << initY << ", initZ = " << initZ << std::endl;
#ifdef SCATTERING
        std::cerr << "scatteredX = " << scatteredX << ", scatteredY = " << scatteredY << ", scatteredZ = " << scatteredZ << std::endl;
#endif
        std::cerr << "hitX = " << hitX << ", hitY = " << hitY << ", hitZ = " << hitZ << std::endl << std::endl; 
      }
    }
#endif

    SVector3 x1(hitX,hitY,hitZ);
    SMatrixSym33 covXYZ = ROOT::Math::SMatrixIdentity();

#ifdef SOLID_SMEAR
    covXYZ(0,0) = varXY; // yn^2 / (xn^2 + yn^2) * delx^2 + xn^2 / (xn^2 + yn^2) * dely^2
    covXYZ(1,1) = varXY; // xn^2 / (xn^2 + yn^2) * delx^2 + yn^2 / (xn^2 + yn^2) * dely^2
    // covXYZ(0,1) -> -(xn * yn) / (xn^2 + yn^2) * delx^2 + (xn * yn) / (xn^2 + yn^2) * dely^2 
    // covXYZ(1,0)  = covXYZ(0,1)    
    covXYZ(2,2) = varZ;
#else
    covXYZ(0,0) = hitX*hitX*varR/hitRad2 + hitY*hitY*varPhi;
    covXYZ(1,1) = hitX*hitX*varPhi + hitY*hitY*varR/hitRad2;
    covXYZ(2,2) = varZ;
    covXYZ(0,1) = hitX*hitY*(varR/hitRad2 - varPhi);
    covXYZ(1,0) = covXYZ(0,1);

    dprint("initPhi: " << initPhi << " hitPhi: " << hitPhi << " initRad: " << initRad  << " hitRad: " << hitRad << std::endl
        << "initX: " << initX << " hitX: " << hitX << " initY: " << initY << " hitY: " << hitY << " initZ: " << initZ << " hitZ: " << hitZ << std::endl 
        << "cov(0,0): " << covXYZ(0,0) << " cov(1,1): " << covXYZ(1,1) << " varZ: " << varZ << " cov(2,2): " << covXYZ(2,2) << std::endl 
        << "cov(0,1): " << covXYZ(0,1) << " cov(1,0): " << covXYZ(1,0) << std::endl);
#endif

    SVector3 initVecXYZ(initX,initY,initZ);
    int index = 10*itrack + ihit;
    Hit initHitXYZ(initVecXYZ,covXYZ,0.,0.,index); 
    MCHitInfo initHitInfo(itrack,simLayer,layer_counts[simLayer]);
    initHits.push_back(initHitXYZ);
    initialhitinfo.push_back(initHitInfo);

    Hit hit1(x1,covXYZ,0.,0.,index);
    hits.push_back(hit1);
    tmpState = propState;

    dprint("initHitId: " //<< initHitXYZ.hitID() << " hit1Id: " << hit1.hitID() <<std::endl
                         << "ihit: " << ihit << " layer: " << simLayer << " counts: " << layer_counts[simLayer]);

    ++layer_counts[simLayer]; // count the number of times passed into layer

  } // end loop over nHitsPerTrack
}


#include <string>
#include <fstream>
#include <sstream>

void setupTrackFromTextFile(SVector3& pos, SVector3& mom, SMatrixSym66& covtrk, HitVec& hits, unsigned int itrack,
			    int& charge, float pt, const Geometry& geom, HitVec& initHits, MCHitInfoVec& initialhitinfo)
{

  //fixme: check also event count

  const float hitposerrXY = 0.01;//assume 100mum uncertainty in xy coordinate
  const float hitposerrZ = 0.1;//assume 1mm uncertainty in z coordinate
  const float hitposerrR = hitposerrXY/10.;

  const float varXY  = hitposerrXY*hitposerrXY;
  const float varZ   = hitposerrZ*hitposerrZ;

  unsigned int nLayers = geom.CountLayers();
  unsigned int layer_counts[nLayers];
  for (unsigned int ilayer=0;ilayer<nLayers;++ilayer){
    layer_counts[ilayer]=0;
  }
  unsigned int nTotHit = nLayers; // can tune this number!
  // to include loopers, and would rather add a break on the code if layer ten exceeded
  // if block BREAK if hit.Layer == theGeom->CountLayers() 
  // else --> if (NMAX TO LOOPER (maybe same as 10?) {break;} else {continue;}
  unsigned int simLayer = 0;
  hits.reserve(nTotHit);
  initHits.reserve(nTotHit);

  // std::ifstream infile("cmssw.simtracks.SingleMu10GeV.10k.eta06z5.txt");
  std::ifstream infile("cmssw.simtracks.SingleMu1GeV.1k.eta06z5.txt");
  // std::ifstream infile("cmssw.simtracks.SingleMu1GeVNoMaterial.1k.eta06z5.txt");
  // std::ifstream infile("cmssw.simtracks.SingleMu06GeV.1k.eta06z5.txt");
  std::string line;
  int countTracks = -1;
  unsigned int countHits   = 0;
  bool gotTrack = false;
  while (std::getline(infile, line)) {

    std::istringstream iss(line);
    std::string type;

    //std::cout << line << std::endl;

    iss >> type;

    if (type=="simTrack") {
      countTracks++;
      if (itrack!=countTracks) continue;
      //std::cout << "countTracks=" << countTracks << std::endl;
      //if (countTracks!=1) continue;
      gotTrack = true;
      float x,y,z,px,py,pz;
      int q;
      iss >> x >> y >> z >> px >> py >> pz >> q;

      pos=SVector3(x,y,z);
      charge = q;
      mom=SVector3(px,py,pz);
      covtrk=ROOT::Math::SMatrixIdentity();
      //initial covariance can be tricky
      for (unsigned int r=0; r<6; ++r) {
	for (unsigned int c=0; c<6; ++c) {
	  if (r==c) {
	    if (r<3) covtrk(r,c)=pow(1.0*pos[r],2);//100% uncertainty on position
	    else covtrk(r,c)=pow(1.0*mom[r-3],2);  //100% uncertainty on momentum
	  } else {
	    covtrk(r,c)=0.;                   //no covariance
	  }
	}
      }

    }

    if (type=="simHit" && gotTrack) {

      countHits++;

      float initX,initY,initZ;
      float r,eta;
      float radl,xi;
      iss >> initX >> initY >> initZ >> r >> eta >> radl >> xi;

      float initPhi = atan2(initY,initX);
      float initRad = sqrt(initX*initX+initY*initY);

      UVector3 init_point(initX,initY,initZ);
      simLayer = geom.LayerIndex(init_point);
      simLayer = 1;//fixme

      float hitZ    = initZ;
      float hitPhi  = initPhi;
      float hitRad  = initRad;
      hitZ    += hitposerrZ*g_gaus(g_gen);
      hitPhi  += ((hitposerrXY/initRad)*g_gaus(g_gen));
      hitRad  += (hitposerrR)*g_gaus(g_gen);

      float hitRad2 = hitRad*hitRad;
      float hitX    = hitRad*cos(hitPhi);
      float hitY    = hitRad*sin(hitPhi);
      
      float varPhi = varXY/hitRad2;
      float varR   = hitposerrR*hitposerrR;

      SVector3 x1(hitX,hitY,hitZ);
      SMatrixSym33 covXYZ = ROOT::Math::SMatrixIdentity();
      covXYZ(0,0) = hitX*hitX*varR/hitRad2 + hitY*hitY*varPhi;
      covXYZ(1,1) = hitX*hitX*varPhi + hitY*hitY*varR/hitRad2;
      covXYZ(0,1) = hitX*hitY*(varR/hitRad2 - varPhi);
      covXYZ(1,0) = covXYZ(0,1);
      covXYZ(2,2) = varZ;
    
      SVector3 initVecXYZ(initX,initY,initZ);
      int index = 10*itrack + countHits;
      Hit initHitXYZ(initVecXYZ,covXYZ,radl,xi,index); 
      MCHitInfo initHitInfo(itrack,simLayer,layer_counts[simLayer]);
      initHits.push_back(initHitXYZ);
      initialhitinfo.push_back(initHitInfo);
      
      Hit hit1(x1,covXYZ,radl,xi,index);
      hits.push_back(hit1);

      ++layer_counts[simLayer];

      if (countHits>=nLayers) return;

    }
    
    // process pair (a,b)
  }

}
