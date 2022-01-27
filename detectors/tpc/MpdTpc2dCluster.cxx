
#include "MpdTpc2dCluster.h"
#include "MpdTpcSector.h"
#include <iostream>
#include <iomanip>
#include <set>

ClassImp(MpdTpc2dCluster);

using namespace std;

//......................................................................

MpdTpc2dCluster::MpdTpc2dCluster(Int_t row, Int_t sec)
   : fFlag(0), fId(0), fSecList(0), fRowList(0), fColList(0), fBktList(0), fAdcList(0), fSector(sec), fADCSum(0),
     fCorrel(0)
{
   fMinBkt = 1000;
   fMaxBkt = -1;
   fMinCol = 1000;
   fMaxCol = -1;
}

//......................................................................

MpdTpc2dCluster::MpdTpc2dCluster(const MpdTpc2dCluster &cl) : TObject((const TObject &)cl)
{
   fFlag   = cl.fFlag;
   fId     = cl.fId;
   fMinBkt = cl.fMinBkt;
   fMaxBkt = cl.fMaxBkt;
   fMinCol = cl.fMinCol;
   fMaxCol = cl.fMaxCol;
   fADCSum = cl.fADCSum;
   fAvgCol = cl.fAvgCol;
   fSigCol = cl.fSigCol;
   fAvgBkt = cl.fAvgBkt;
   fSigBkt = cl.fSigBkt;
   fCorrel = cl.fCorrel;

   for (UInt_t i = 0; i < cl.GetNumDigits(); ++i) {
      fSecList.push_back(cl.fSecList[i]);
      fRowList.push_back(cl.fRowList[i]);
      fColList.push_back(cl.fColList[i]);
      fBktList.push_back(cl.fBktList[i]);
      fAdcList.push_back(cl.fAdcList[i]);
   }
}

//......................................................................

MpdTpc2dCluster::~MpdTpc2dCluster() {}

//......................................................................

Bool_t MpdTpc2dCluster::Insert(Int_t row, Int_t col, Int_t bkt, Float_t adc)
{
   Int_t   thiscol, thisbkt;
   Float_t adcval;

   thiscol = col;
   thisbkt = bkt;
   adcval  = adc;

   fRowList.push_back(row);
   fColList.push_back(col);
   fBktList.push_back(bkt);
   fAdcList.push_back(adc);

   // now set the other variables, if necessary

   if (thiscol < fMinCol) fMinCol = thiscol;
   if (thiscol > fMaxCol) fMaxCol = thiscol;
   if (thisbkt < fMinBkt) fMinBkt = thisbkt;
   if (thisbkt > fMaxBkt) fMaxBkt = thisbkt;

   fADCSum += adcval;

   return kTRUE;
}

//......................................................................

Int_t MpdTpc2dCluster::Row(int i) const
{
   return fRowList[i];
}
Int_t MpdTpc2dCluster::Col(int i) const
{
   return fColList[i];
}
Int_t MpdTpc2dCluster::Bkt(int i) const
{
   return fBktList[i];
}
Float_t MpdTpc2dCluster::Adc(int i) const
{
   return fAdcList[i];
}
Int_t MpdTpc2dCluster::Sec(int i) const
{
   return fSecList[i];
}

//......................................................................

Bool_t MpdTpc2dCluster::Insert(Int_t sec, Int_t row, Int_t col, Int_t bkt, Float_t adc)
{
   // Insert pixel

   Bool_t asd = Insert(row, col, bkt, adc);
   fSecList.push_back(sec);
   return asd;
}

//......................................................................

Int_t MpdTpc2dCluster::NTracks() const
{
   // Get number of track contributors to this cluster

   std::set<Int_t> ids;
   Int_t           ndig = fSecList.size();
   for (Int_t i = 0; i < ndig; ++i) ids.insert(fSecList[i]);
   return ids.size();
}

//......................................................................

////////////////////////////////////////////////////////////////////////