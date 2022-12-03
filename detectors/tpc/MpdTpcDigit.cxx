#include "MpdTpcDigit.h"
using std::map;

ClassImp(MpdTpcDigit);

MpdTpcDigit::MpdTpcDigit() : BaseTpcDigit()
{
   fOrigin = -1;
}

MpdTpcDigit::MpdTpcDigit(Int_t ori, Int_t pad, Int_t row, Int_t bin, Int_t sec, Float_t adc)
   : BaseTpcDigit(sec, row, pad, bin, adc)
{
   fOrigin = ori;
}

MpdTpcDigit::~MpdTpcDigit() {}
