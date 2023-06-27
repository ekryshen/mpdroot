//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_FILE_SOURCE_H
#define __MPD_FILE_SOURCE_H

//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFileSource
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------

#define private protected
#include "FairFileSource.h"
#undef private
//------------------------------------------------------------------------------------------------------------------------
class MpdFileSource : public FairFileSource 
{

public:
    	MpdFileSource(const char* flnm);

	virtual Bool_t 	Init();
	virtual void 	FillEventHeader(FairEventHeader*);

    ClassDef(MpdFileSource, 1)
};
//------------------------------------------------------------------------------------------------------------------------
#endif 
