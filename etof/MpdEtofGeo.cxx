//------------------------------------------------------------------------------------------------------------------------
#include "FairGeoNode.h"

#include "MpdEtofGeo.h"
//------------------------------------------------------------------------------------------------------------------------
MpdEtofGeo::MpdEtofGeo() 
{
	fName="etof";
	maxSectors=0;
	maxModules=24;
}
//------------------------------------------------------------------------------------------------------------------------
const char* MpdEtofGeo::getModuleName(Int_t m) 
{
	sprintf(modName,"etof%i",m+1);
return modName;
}
//------------------------------------------------------------------------------------------------------------------------
const char* MpdEtofGeo::getEleName(Int_t m) 
{
	sprintf(eleName,"et%i",m+1);
return eleName;
}
//------------------------------------------------------------------------------------------------------------------------
ClassImp(MpdEtofGeo)
