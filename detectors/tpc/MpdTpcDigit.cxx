#include "MpdTpcDigit.h"
using std::map;

ClassImp(MpdTpcDigit);

MpdTpcDigit::MpdTpcDigit() : fPad(-1), fRow(-1), fTimeBin(-1), fSector(-1), fAdc(-1.0), fOrigin(-1), fOrigWeight(-1) {}

//......................................................................

MpdTpcDigit::MpdTpcDigit(Int_t ori, Int_t pad, Int_t row, Int_t bin, Int_t sec, Float_t adc)
   : fPad(pad), fRow(row), fTimeBin(bin), fSector(sec), fAdc(adc), fOrigin(ori), fOrigWeight(-1)
{
}

//......................................................................

MpdTpcDigit::MpdTpcDigit(Int_t ori, Int_t pad, Int_t row, Int_t bin, Int_t sec, Float_t adc,
                         std::map<Int_t, Float_t> &origins)
   : fPad(pad), fRow(row), fTimeBin(bin), fSector(sec), fAdc(adc), fOrigin(ori), fOrigWeight(-1)
{
   fOrigWeight = origins[ori];
   fOrigins.insert(origins.begin(), origins.end());
}

//......................................................................
/*AZ
MpdTpcDigit::MpdTpcDigit(MpdTpcDigit* digit) :
fRow(digit->GetRow()), fPad(digit->GetPad()), fOrigin(digit->GetOrigin()), fSector(digit->GetSector()),
fAdc(digit->GetAdc()), fTimeBin(digit->GetTimeBin()) {}
*/
//......................................................................

MpdTpcDigit::~MpdTpcDigit() {}
