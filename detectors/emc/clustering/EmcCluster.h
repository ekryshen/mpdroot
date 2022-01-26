// $Id$
// Author: artur   2016/04/15

#ifndef __EMCCLUSTER_H__
#define __EMCCLUSTER_H__

#include "Cluster.h"
#include "EmcClusterElement.h"
#include "EmcClusterInfo.h"

class TH2D;

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// EmcCluster                                                                 //
//                                                                            //
// <brief class description>                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

class EmcCluster: public Cluster {

public:

    EmcCluster();
    EmcCluster(Int_t uid);
    
    virtual ~EmcCluster();
    
    EmcClusterInfo* CreateClusterInfo();
    
    virtual void FillClusterInfo(EmcClusterInfo* clinfo);
    virtual void FillH(TH2D* h);

private:

    ClassDef(EmcCluster,1)
};

#endif  /* __EMCCLUSTER_H__ */

