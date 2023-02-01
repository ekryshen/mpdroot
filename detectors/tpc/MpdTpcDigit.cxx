#include "MpdTpcDigit.h"

ClassImp(MpdTpcDigit);

MpdTpcDigit::MpdTpcDigit() : fPad(-1), fRow(-1), fTimeBin(-1), fSector(-1), fAdc(-1.0), fOrigin(-1) {}

//......................................................................

MpdTpcDigit::MpdTpcDigit(Int_t ori, Int_t pad, Int_t row, Int_t bin, Int_t sec, Float_t adc)
   : fPad(pad), fRow(row), fTimeBin(bin), fSector(sec), fAdc(adc), fOrigin(ori)
{
}

MpdTpcDigit::~MpdTpcDigit() {}
