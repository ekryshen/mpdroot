//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_TOF_POINT_H
#define __MPD_TOF_POINT_H 1

#include <assert.h>

#include "TObject.h"
#include "TVector3.h"

#include "FairMCPoint.h"
//------------------------------------------------------------------------------------------------------------------------
class MpdTofPoint : public FairMCPoint {
public:
   /** Constructor with arguments
    *@param tid      Index of MCTrack
    *@param suid     strip unique ID
    *@param pos      Ccoordinates at entrance to active volume [cm]
    *@param mom      Momentum of track at entrance [GeV]
    *@param tof      Time since event start [ns]
    *@param length   Track length since creation [cm]
    *@param eLoss    Energy deposit [GeV]
    **/
   MpdTofPoint(Int_t tid, Int_t suid, TVector3 pos, TVector3 mom, Double_t tof, Double_t length, Double_t eLoss);
   MpdTofPoint();
   virtual ~MpdTofPoint(){};

   // Output to screen
   virtual void Print(const Option_t *opt) const;

   // CAUTION:  MAX_VALUE = 255(0xFF)

   // MpdTof --------------------------------
   // sector 	[1,...,14],	0xFF000000
   // detector	[1,...,20],	0x00FF0000
   // gap	 	[1,...,3],	0x0000FF00
   // strip 	[1,...,24],	0x000000FF

   // MpdEtof --------------------------------
   // strip 	[1,...,131],??	0x000000FF
   // detector 	[1,...,1],??	0x0000FF00
   // box		[1,...,30],??	0x00FF0000
   // sector 	[1,...,2],??	0xFF000000

   inline Int_t GetSector() const { return ((fDetectorID & 0xFF000000) >> 24); };
   inline Int_t GetDetector() const { return ((fDetectorID & 0x00FF0000) >> 16); };
   inline Int_t GetGap() const { return ((fDetectorID & 0x0000FF00) >> 8); };
   inline Int_t GetStrip() const { return (fDetectorID & 0x000000FF); };

   static Int_t GetSector(Int_t suid) { return ((suid & 0xFF000000) >> 24); };
   static Int_t GetDetector(Int_t suid) { return ((suid & 0x00FF0000) >> 16); };
   static Int_t GetGap(Int_t suid) { return ((suid & 0x0000FF00) >> 8); };
   static Int_t GetStrip(Int_t suid) { return (suid & 0x000000FF); };

   static Int_t ClearGap(Int_t suid) { return (suid & 0xFFFF00FF); };                // set gap = 0 for hit suid
   static Int_t SetCentralGap(Int_t suid) { return (ClearGap(suid) | 0x00000200); }; // set gap = 2 for hit suid

   Int_t       GetSuid() const { return fDetectorID; };
   static bool PrintSuid(Int_t suid, const char *comment = nullptr, std::ostream &os = std::cout);

   static bool  IsSameDetector(Int_t suid1, Int_t suid2);
   static bool  IsSameStrip(Int_t suid1, Int_t suid2);
   static bool  IsSameGap(Int_t suid1, Int_t suid2);
   static void  ParseSuid(Int_t suid, Int_t &sector, Int_t &detector, Int_t &gap, Int_t &strip);
   static Int_t GetSuid72(Int_t sector, Int_t detector, Int_t strip);            // stripID [1,72]
   static Int_t GetSuid24(Int_t sector, Int_t detector, Int_t gap, Int_t strip); // stripID [1,24]

   ClassDef(MpdTofPoint, 2);
};
//------------------------------------------------------------------------------------------------------------------------
#endif
