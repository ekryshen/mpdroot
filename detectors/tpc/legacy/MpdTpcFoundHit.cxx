#include "MpdTpcFoundHit.h"

ClassImp(MpdTpcFoundHit);

MpdTpcFoundHit::MpdTpcFoundHit()
   : fNumHits(1), fHitType(kUnknown), fPadRow(0), fPadCol(0), fTimeBkt(0), fQFit(0.), fSigQFit(0.), fQADC(0.),
     fCluster(0)
{
}

//......................................................................

MpdTpcFoundHit::MpdTpcFoundHit(const MpdTpc2dCluster *clus)
   : fNumHits(1), fHitType(kUnknown), fPadRow(0), fPadCol(0), fTimeBkt(0), fQFit(0.), fSigQFit(0.), fQADC(0.),
     fCluster(0)
{
   this->AttachCluster(clus);
}

//......................................................................

MpdTpcFoundHit::MpdTpcFoundHit(Float_t x, Float_t dx, Float_t y, Float_t dy, Float_t z, Float_t dz)
{
   fNumHits = 1;
   fHitType = kUnknown;
   fQFit    = 0.;
   fSigQFit = 0.;
   fQADC    = 0.;
   fTimeBkt = 0.;

   Float_t pos[3] = {x, y, z};
   Float_t err[3] = {dx, dy, dz};
   this->SetPos(pos, err);
}

//......................................................................

MpdTpcFoundHit::MpdTpcFoundHit(Float_t p[3], Float_t dp[3])
{
   fNumHits = 1;
   fHitType = kUnknown;
   fQFit    = 0.;
   fSigQFit = 0.;
   fQADC    = 0.;
   fTimeBkt = 0.;
   this->SetPos(p, dp);
}

//......................................................................

MpdTpcFoundHit::MpdTpcFoundHit(const MpdTpcFoundHit &hit) : TObject((const TObject &)hit)
{

   this->AttachCluster(hit.Cluster());
   fTimeBkt = hit.fTimeBkt;
   fNumHits = hit.fNumHits;
   fHitType = hit.fHitType;
   fQFit    = hit.fQFit;
   fSigQFit = hit.fSigQFit;
   fQADC    = hit.fQADC;
   fPadCol  = hit.fPadCol;
}

//......................................................................

MpdTpcFoundHit::~MpdTpcFoundHit() {}

const void MpdTpcFoundHit::Print() const
{
   std::cout.precision(4);
   std::cout << "(" << this->GetLocalX() << " +/- " << this->errX() << "," << this->GetLocalY() << " +/- "
             << this->errY() << "," << this->GetLocalZ() << " +/- " << this->errZ() << ")";
}

//......................................................................

void MpdTpcFoundHit::AttachCluster(const MpdTpc2dCluster *currclus)
{
   // fCluster = TRef((MpdTpc2dCluster*)currclus);
   fCluster = currclus;
   fPadRow  = currclus->Row();
}

//......................................................................

ostream &operator<<(ostream &ostr, const MpdTpcFoundHit &h1)
{
   ostr.precision(4);
   ostr << "(" << h1.GetLocalX() << "," << h1.GetLocalY() << "," << h1.GetLocalZ() << ")";
   return ostr;
}

void MpdTpcFoundHit::SetPos(Float_t *p, Float_t *dp)
{
   _xl = p[0];
   _yl = p[1];
   _zl = p[2];
   _dx = dp[0];
   _dy = dp[1];
   _dz = dp[2];
}

void MpdTpcFoundHit::AddOrigin(Int_t Origin)
{
   //    for(Int_t i = 0; i < Origins.GetEntriesFast(); ++i){
   //	Bool_t OriginExist = 0;
   //	for(Int_t j = 0; j < fOrigins.GetEntriesFast(); ++j)
   //	    if(fOrigins.At(j) == Origins.At(i)) {
   //                OriginExist = 1;
   //                continue;
   //            }
   //	if(!OriginExist)
   //	    fOrigins.Add(Origins[i]);
   //    }
   fOrigin = Origin;
}

////////////////////////////////////////////////////////////////////////
