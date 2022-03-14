#include "TpcLheHit.h"
#include "TpcLheTrack.h"
#include "TClonesArray.h"

#include "Riostream.h"
#include <iomanip>

using namespace std;

////////////////////////////////////////////////////////////////////////
//                                                                    //
// TpcLheTrack  class - track representation for the LHE              //
//                                                                    //
////////////////////////////////////////////////////////////////////////

ClassImp(TpcLheTrack)

   //______________________________________________________________
   TpcLheTrack::TpcLheTrack()
{
   // Default constructor.

   SetDefaults();
}

//______________________________________________________________
TpcLheTrack::TpcLheTrack(Int_t tracknumber)
{
   // the track number is set.

   //  cout << " TpcLheTrack (tracknumber)" << endl;

   SetDefaults();
   SetTrackNumber(tracknumber);
}

//______________________________________________________________
TpcLheTrack ::~TpcLheTrack()
{
   // Destructor.

   //  cout << " Destructor for LheTrack" << endl;
   if (fRealHits) { // fRealHits->Delete();} //;
      delete fRealHits;
      fRealHits = NULL;
   }
}

//___________________________________________________________
void TpcLheTrack ::SetDefaults()
{
   // Default setup for the track.

   //  fRealHits = new TRefArray();
   fRealHits = new TObjArray(0);

   fTrackNumber = -1;

   SetPid(0);

   ComesFromMainVertex(kFALSE);

   fP.SetX(0.);
   fP.SetY(0.);
   fP.SetZ(0.);

   fVertex.SetX(0.);
   fVertex.SetY(0.);
   fVertex.SetZ(0.);

   fQ = 0;
}

//______________________________________________________________
void TpcLheTrack ::SetTrackNumber(Int_t number)
{

   // Sets the tracknumber for track.
   fTrackNumber = number;
}

//______________________________________________________________
void TpcLheTrack::SetVertex(Double_t vx, Double_t vy, Double_t vz)
{

   fVertex.SetX(vx);
   fVertex.SetY(vy);
   fVertex.SetZ(vz);
}

//______________________________________________________________
void TpcLheTrack::SetCircle(Double_t x, Double_t y, Double_t r)
{

   fCircle.SetX(x);
   fCircle.SetY(y);
   fCircle.SetZ(r);
}

//______________________________________________________________
void TpcLheTrack::AddHit(TpcLheHit *point)
{
   // Adds a hit to the track.

   fRealHits->AddLast(point);
   //  point->SetUsage(kTRUE);
}

//______________________________________________________________
void TpcLheTrack::PrintHits()
{
   //

   Int_t nhit = fRealHits->GetEntries();

   TpcLheHit *ghit = NULL;

   for (Int_t j = 0; j < nhit; j++) {
      ghit = (TpcLheHit *)fRealHits->At(j);
      cout << " " << ghit->GetTrackID();
   }
   cout << endl;
}

//______________________________________________________________
void TpcLheTrack::Print()
{
   //

   Int_t nhit = fRealHits->GetEntries();

   cout << "\n Track " << fTrackNumber << "  nhits " << nhit;
   cout << " pdg " << fPid;
   cout << " primary " << fFromMainVertex << endl;

   cout << " vertex XYZ: " << setw(5) << fVertex.GetX();
   cout << " " << setw(5) << fVertex.GetY();
   cout << " " << setw(5) << fVertex.GetZ();

   cout << " Pxyz: ";
   cout.precision(3);
   cout << " " << setw(7) << fP.X();
   cout << " " << setw(7) << fP.Y();
   cout << " " << setw(7) << fP.Z();
   cout << endl;

   TpcLheHit *ghit = NULL;

   for (Int_t j = 0; j < nhit; j++) {
      ghit = (TpcLheHit *)fRealHits->At(j);
      ghit->Print();
   }
}
