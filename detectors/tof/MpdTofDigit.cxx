//--------------------------------------------------------------------------------------------------------------------------------------
#include "MpdTofDigit.h"

ClassImp(MpdTofDigit)

//--------------------------------------------------------------------------------------------------------------------------------------
MpdTofDigit::MpdTofDigit(Short_t sector, Short_t plane, Short_t strip, Short_t side,Float_t t,Float_t a)
: fSector(sector), fPlane(plane), fStrip(strip), fSide(side), fTime(t), fAmplitude(a)
{

}
//--------------------------------------------------------------------------------------------------------------------------------------
MpdTofDigit::MpdTofDigit(const MpdTofDigit *ptr, Float_t t, Float_t a)
: fPlane(ptr->fPlane), fStrip(ptr->fStrip), fSide(ptr->fSide), fTime(t), fAmplitude(a)
{

}
//--------------------------------------------------------------------------------------------------------------------------------------
void	MpdTofDigit::print(const char* comment, std::ostream& os)const
{
	if(nullptr != comment) os<<comment;

	os<<"  sector: " << fSector << ", det: "<<fPlane<<", strip: "<<fStrip<<", stripSide: "<<fSide<<", time: "<<fTime<<", ampl.: "<<fAmplitude;
}
//--------------------------------------------------------------------------------------------------------------------------------------

