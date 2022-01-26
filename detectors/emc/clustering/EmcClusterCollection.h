// $Id$
// Author: artur   2016/04/15

#ifndef __EMCCLUSTERCOLLECTION_H__
#define __EMCCLUSTERCOLLECTION_H__

#include <TObject.h>
#include "ClusterCollection.h"
#include "EmcCluster.h"

class TClonesArray;
class TH2D;

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// EmcClusterCollection                                                       //
//                                                                            //
// <brief class description>                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

class EmcClusterCollection: public ClusterCollection {

public:

    EmcClusterCollection();
    EmcClusterCollection(Int_t uid);
    //EmcClusterCollection(TString name);
    
    virtual ~EmcClusterCollection();

    virtual Bool_t FillClusterInfo(TClonesArray* info, Int_t clflag = -1);
    virtual void   FillH(TH2D* h);
       
private:

    ClassDef(EmcClusterCollection,1)
};

#endif  /* __EMCCLUSTERCOLLECTION_H__ */

