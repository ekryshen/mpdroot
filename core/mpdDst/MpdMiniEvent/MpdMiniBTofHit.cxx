//
// MpdMiniBTofHit holds information about BTof hits
//

// MiniDst headers
#include "MpdMiniMessMgr.h"
#include "MpdMiniBTofHit.h"

ClassImp(MpdMiniBTofHit)

   //_________________
   MpdMiniBTofHit::MpdMiniBTofHit()
   : TObject(), fDetectorID(-1), fBTofHitPosX(0), fBTofHitPosY(0), fBTofHitPosZ(0), fTime(0)
{
   // Default constructor
}

//_________________
MpdMiniBTofHit::MpdMiniBTofHit(int id) : MpdMiniBTofHit()
{
   // Parametrized constructor
   if (id < 0) return;
   fDetectorID = ((id > std::numeric_limits<short>::max()) ? std::numeric_limits<short>::max() : (Short_t)id);
}

//_________________
MpdMiniBTofHit::MpdMiniBTofHit(const MpdMiniBTofHit &hit) : TObject()
{
   // Copy constructor
   fDetectorID  = hit.fDetectorID;
   fTime        = hit.fTime;
   fBTofHitPosX = hit.fBTofHitPosX;
   fBTofHitPosY = hit.fBTofHitPosY;
   fBTofHitPosZ = hit.fBTofHitPosZ;
}

//_________________
MpdMiniBTofHit::~MpdMiniBTofHit()
{
   /* empty */
}

//_________________
void MpdMiniBTofHit::Print(const Char_t *option __attribute__((unused))) const
{
   LOG_INFO << " Id: " << id() << " x/y/z: " << btofHitPosX() << "/" << btofHitPosY() << "/" << btofHitPosZ()
            << " time: " << time() << endm;
}

//_________________
void MpdMiniBTofHit::setId(Int_t sector, Int_t gap, Int_t detector, Int_t strip)
{
   fDetectorID = (sector << 24) + (detector << 16) + (gap << 8) + strip;
}

//_________________
void MpdMiniBTofHit::setHitPositionX(Float_t x)
{
   fBTofHitPosX = (fabs(x * 100.) > std::numeric_limits<short>::max()
                      ? ((x > 0) ? std::numeric_limits<short>::max() : std::numeric_limits<short>::min())
                      : (Short_t)(TMath::Nint(x * 100.)));
}

//_________________
void MpdMiniBTofHit::setHitPositionY(Float_t y)
{
   fBTofHitPosY = (fabs(y * 100.) > std::numeric_limits<short>::max()
                      ? ((y > 0) ? std::numeric_limits<short>::max() : std::numeric_limits<short>::min())
                      : (Short_t)(TMath::Nint(y * 100.)));
}

//_________________
void MpdMiniBTofHit::setHitPositionZ(Float_t z)
{
   fBTofHitPosZ = (fabs(z * 100.) > std::numeric_limits<short>::max()
                      ? ((z > 0) ? std::numeric_limits<short>::max() : std::numeric_limits<short>::min())
                      : (Short_t)(TMath::Nint(z * 100.)));
}

//_________________
void MpdMiniBTofHit::setHitPositionXYZ(Float_t x, Float_t y, Float_t z)
{
   setHitPositionX(x);
   setHitPositionY(y);
   setHitPositionZ(z);
}
