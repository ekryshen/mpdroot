//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFwdTrackProducer
/// \brief Ideal track producer
/// \author Evgeny Kryshen (PNPI, Gatchina)
//------------------------------------------------------------------------------------------------------------------------
#include "map"
#include "unordered_set"

#include "TClonesArray.h"
#include "TF1.h"
#include "TGraph.h"

#include "FairRootManager.h"
#include "MpdFwdHit.h"
#include "MpdFwdTrack.h"
#include "MpdFwdPoint.h"
#include "MpdMCTrack.h"
#include "MpdFwdTrackProducer.h"

ClassImp(MpdFwdTrackProducer)

MpdFwdTrackProducer::MpdFwdTrackProducer() : FairTask("FwdTrackProducer"){}
MpdFwdTrackProducer::~MpdFwdTrackProducer(){}

InitStatus MpdFwdTrackProducer::Init(){
  fMcTracks = (TClonesArray*) FairRootManager::Instance()->GetObject("MCTrack");
  fFwdPoints = (TClonesArray*) FairRootManager::Instance()->GetObject("FwdPoint");
  fFwdHits = (TClonesArray*) FairRootManager::Instance()->GetObject("FwdHit");
  fFwdTracks = new TClonesArray("MpdFwdTrack");
  FairRootManager::Instance()->Register("FwdTrack", "Fwd", fFwdTracks, kTRUE);
  return kSUCCESS;
}

void MpdFwdTrackProducer::Exec(Option_t *opt){
  FindTracks();
  FitTracks();
  PrintTracks();
}


void MpdFwdTrackProducer::FindTracks(){ 
  // clear reconstructed track array
  fFwdTracks->Clear();
  
  // map from MC track to reconstructed track
  std::unordered_map<Int_t, Int_t> mapFwdTracks;
  
  // Collect hits using MC info and create corresponding reconstructed tracks
  Int_t itrk = 0;
  for (Int_t ihit=0; ihit<fFwdHits->GetEntriesFast(); ihit++){
    MpdFwdHit* hit = (MpdFwdHit*) fFwdHits->UncheckedAt(ihit);
    MpdFwdPoint* point = (MpdFwdPoint*) fFwdPoints->UncheckedAt(hit->GetRefIndex());
    Int_t trackID = point->GetTrackID();
    MpdMCTrack* mcTrack = (MpdMCTrack*) fMcTracks->UncheckedAt(trackID);
    if (mcTrack->GetMotherId()>=0) continue; // skip secondary tracks for the moment
    MpdFwdTrack* tr = nullptr;
    if (mapFwdTracks.find(trackID)==mapFwdTracks.end()){
      // no hits for this track yet -> create new
      mapFwdTracks[trackID] = itrk;
      tr = new ((*fFwdTracks)[itrk++]) MpdFwdTrack();
    } else {
      // track is found
      tr = (MpdFwdTrack*) fFwdTracks->UncheckedAt(mapFwdTracks.at(trackID)); 
    }
    tr->AddHitIndex(ihit);
  }
}

void MpdFwdTrackProducer::FitTracks(){
  TF1 funcCircle("funcCircle","[3]*sqrt([0]^2-(x-[1])^2)+[2]",-100,100);
  TGraph g; // graph of X,Y coordinates
  const double B  = 0.5; // T
  int q;
  double x1,x2,x3,y1,y2,y3,ma,mb,r0,x0,y0,r,xc,yc,x,y,tx,ty,tz,xl,yl,zl,pt,phi,lr;
  double v[6]; // state vector
  
  for (Int_t itrack=0;itrack<fFwdTracks->GetEntriesFast();itrack++){
    MpdFwdTrack* track = (MpdFwdTrack*) fFwdTracks->UncheckedAt(itrack);
    Int_t n= track->GetNIndices();
    if (n<3) continue; // unable to fit tracks with less than 3 hits
    g.Set(n);
    for (Int_t ihit = 0; ihit < n; ihit++){
      Int_t hitIndex = track->GetHitIndex(ihit);
      MpdFwdHit* hit = (MpdFwdHit*) fFwdHits->UncheckedAt(hitIndex);
      g.GetX()[ihit] = hit->GetX();
      g.GetY()[ihit] = hit->GetY();
      if (ihit != n-1) continue;
      // last measured hit
      xl = hit->GetX();
      zl = hit->GetZ();
    }
    
    // deduce initial circle parameters using first, middle and last point 
    x1 = g.GetX()[  0]; 
    y1 = g.GetY()[  0]; 
    x2 = g.GetX()[n/2]; 
    y2 = g.GetY()[n/2]; 
    x3 = g.GetX()[n-1]; 
    y3 = g.GetY()[n-1]; 
    ma = (y2-y1)/(x2-x1);
    mb = (y3-y2)/(x3-x2);
    x0 = (ma*mb*(y1-y3)+mb*(x1+x2)-ma*(x2+x3))/2./(mb-ma); 
    y0 = (y1+y2)/2.-(x0-(x1+x2)/2.)/ma;
    r0 = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
    funcCircle.SetParameters(r0,x0,y0,1.);
    
    // fix sign (orientation of the circle related to track charge)
    funcCircle.FixParameter(3,1.);
    if (fabs(funcCircle.Eval(x1)-y1)>1e-3) funcCircle.FixParameter(3,-1);
    if (fabs(funcCircle.Eval(x1)-y1)>1e-3) {
      printf("ERROR\n");
      continue;
    }
    
    // fit track using circle equation
    g.Fit(&funcCircle,"Q");
    r  = funcCircle.GetParameter(0);
    xc = funcCircle.GetParameter(1);
    yc = funcCircle.GetParameter(2);
    q  = round(funcCircle.GetParameter(3));
    
    // track parameters at z=0
    x  = 0;
    y  = funcCircle.Eval(x); // y coordinate at x=0. 
    ty = -(x-xc)/(y-yc); // ty = py/px
    pt = 0.3*B*r/100; // 100 to convert from cm to m
    
    // deduce radial track length from track parameters at last measured z
    yl = funcCircle.Eval(xl);
    phi = acos(1-(xl*xc+yl*yc)/r/r);
    lr = r*phi;
    tz = lr/zl; // tz = pt/pz
    
    // fill state vector (z, x, y, q/pt, ty, tz)
    v[0] = 0;
    v[1] = x;
    v[2] = y;
    v[3] = ty;
    v[4] = tz;
    v[5] = q/pt;
    track->SetStateVector(v);
  }
}

void MpdFwdTrackProducer::PrintTracks(){
  int nTracks = fFwdTracks->GetEntriesFast();
  printf("N tracks: %d\n", nTracks);
  for (Int_t itr=0; itr<fFwdTracks->GetEntriesFast(); itr++){
    MpdFwdTrack* tr = (MpdFwdTrack*) fFwdTracks->UncheckedAt(itr);
    tr->Print();
    printf("Hit z-positions: ");
    for (Int_t ihit=0;ihit<tr->GetNIndices();ihit++){
      Int_t hitIndex = tr->GetHitIndex(ihit);
      MpdFwdHit* hit = (MpdFwdHit*) fFwdHits->UncheckedAt(hitIndex);
      printf("%.2f ",hit->GetZ());
    }
    printf("\n");
  }
}

void MpdFwdTrackProducer::Finish(){
}
