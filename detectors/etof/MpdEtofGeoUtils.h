//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_ETOF_GEO_UTILS_H
#define __MPD_ETOF_GEO_UTILS_H 1

#include <TVector3.h>
#include <TList.h>

#include "IntervalTree.h"
#include "MpdTofGeoUtils.h"

//------------------------------------------------------------------------------------------------------------------------
class TH1D;
class TH2D;
class MpdEtofGeoUtils {
   TmStrips      mStrips;          //! indexing strips by detectorUID
   TintervalTree mDetectorsR[2];   //! detectors R[cm] location interval tree, two sides
   TintervalTree mDetectorsPhi[2]; //! detectors Phi[rads] location interval tree, two sides

   MpdEtofGeoUtils();
   static MpdEtofGeoUtils *instance;

public:
   static MpdEtofGeoUtils *Instance()
   {
      if (instance == nullptr) instance = new MpdEtofGeoUtils;
      return instance;
   }

   void          FindNeighborStrips(Double_t thresh,
                                    bool     forced = false); // thresh, [cm] <--- thresh. distance between neighbor strips
   void          ParseTGeoManager(bool forced = false);
   const LStrip *FindStrip(Int_t UID);

   const TintervalTree *GetDetR() { return mDetectorsR; };
   const TintervalTree *GetDetPhi() { return mDetectorsPhi; };
};
//------------------------------------------------------------------------------------------------------------------------
#endif
