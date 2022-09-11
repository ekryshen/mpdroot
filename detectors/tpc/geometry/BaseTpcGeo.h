//-----------------------------------------------------------
// Description:
//      BaseTpcGeo class is the base class for TPC Geometry.
//      Any custom Tpc geometry must be its child class.
//
// Environment:
//      Software for the MPD Detector at NICA.
//
// Author:
//      Slavomir Hnatic LIT, JINR, Dubna - 9.2022
//
//-----------------------------------------------------------

#ifndef BASETPCGEO_HH
#define BASETPCGEO_HH

#include <TObject.h>

// FairRoot Class Headers ----------------

// MpdRoot Class Headers -----------------

// ROOT Class Declarations ---------------

class BaseTpcGeo : public TObject {
public:
   // Constructors/Destructors ---------
   BaseTpcGeo();
   virtual ~BaseTpcGeo();

protected:

private:

   ClassDef(BaseTpcGeo, 1);
};

#endif


