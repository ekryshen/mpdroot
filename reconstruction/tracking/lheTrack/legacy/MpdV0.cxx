#include "MpdV0.h"

// Default constructor
MpdV0::MpdV0()
   : TNamed("V0", "Global"), fID0(0), fID1(0), fID2(0), fB(0), fchi2(0), fmass(0), fcosa(0), fDCA(0), fPx(0), fPy(0),
     fPz(0), fRx(0), fRy(0), fRz(0), fX(0), fY(0), fZ(0)
{
}
//__________________________________________________________________________

// Constructor with name and title

MpdV0::MpdV0(const char *name, const char *title)
   : TNamed(name, title), fID0(0), fID1(0), fID2(0), fB(0), fchi2(0), fmass(0), fcosa(0), fDCA(0), fPx(0), fPy(0),
     fPz(0), fRx(0), fRy(0), fRz(0), fX(0), fY(0), fZ(0)
{
}
//__________________________________________________________________________

// constructor with all parameters
MpdV0::MpdV0(const char *name, const char *title, Double_t vtx[3], Double_t chi2, Double_t p[3], Double_t DCA,
             Double_t cosa, Double_t mass, Double_t r[3], Int_t B, Int_t id1, Int_t id2, Int_t id0)
   : TNamed(name, title)
{
   fID0  = id0;
   fID1  = id1;
   fID2  = id2;
   fB    = B;
   fchi2 = chi2;
   fmass = mass;
   fcosa = cosa;
   fDCA  = DCA;
   fPx   = p[0];
   fPy   = p[1];
   fPz   = p[2];
   fRx   = r[0];
   fRy   = r[1];
   fRz   = r[2];
   fX    = vtx[0];
   fY    = vtx[1];
   fZ    = vtx[2];
}
//__________________________________________________________________________
MpdV0::~MpdV0(){};

void MpdV0::SetV0(Double_t vtx[3], Double_t chi2, Double_t p[3], Double_t DCA, Double_t cosa, Double_t mass,
                  Double_t r[3], Int_t B, Int_t id1, Int_t id2, Int_t id0)
{
   fID0  = id0;
   fID1  = id1;
   fID2  = id2;
   fB    = B;
   fchi2 = chi2;
   fmass = mass;
   fcosa = cosa;
   fDCA  = DCA;
   fPx   = p[0];
   fPy   = p[1];
   fPz   = p[2];
   fRx   = r[0];
   fRy   = r[1];
   fRz   = r[2];
   fX    = vtx[0];
   fY    = vtx[1];
   fZ    = vtx[2];
}

MpdV0 &MpdV0::operator=(const MpdV0 &V0)
{
   if (this == &V0) return *this;
}

Double_t *MpdV0::GetPxPyPz(Double_t p[3])
{
   p[0] = fPx;
   p[1] = fPy;
   p[2] = fPz;
   return p;
}

Double_t *MpdV0::GetXYZ(Double_t vtx[3])
{
   vtx[0] = fX;
   vtx[1] = fY;
   vtx[2] = fZ;
   return vtx;
}

Double_t MpdV0::GetP()
{
   Double_t P = TMath::Sqrt(fPx * fPx + fPy * fPy + fPz * fPz);
   return P;
}

Double_t MpdV0::GetDist()
{
   Double_t dist = TMath::Sqrt(fRx * fRx + fRy * fRy + fRz * fRz);
   return dist;
}

ClassImp(MpdV0);