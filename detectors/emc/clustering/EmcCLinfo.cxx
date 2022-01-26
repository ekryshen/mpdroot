// $Id$
// Author: artur   2016/04/25


//_____________________________________________________________________________
//
// EmcCLinfo
//_____________________________________________________________________________

#include "EmcCLinfo.h"
#include <fstream>
#include <iostream>

using namespace std;

ClassImp(EmcCLinfo)

//_____________________________________________________________________________
EmcCLinfo::EmcCLinfo():fNEvents(0),fCLMethod(0)
{
 
}

//_____________________________________________________________________________
EmcCLinfo::~EmcCLinfo() 
{
 
}

//_____________________________________________________________________________
void EmcCLinfo::Print(const Option_t* opt) const
{
  
}

//_____________________________________________________________________________
void EmcCLinfo::Save(TString fullname)
{
   ofstream ff(fullname.Data());
   if (!ff.good()) return;
   
   ff << "<EmcCLinfo>" << endl;
   ff << "Processed events: " << fNEvents << endl;
   ff << "Clusterization method: " << fCLMethod << endl;
   
   ff.close();
}






