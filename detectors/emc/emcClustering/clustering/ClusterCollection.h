// $Id$
// Author: artur   2016/04/07

#ifndef __CLUSTERCOLLECTION_H__
#define __CLUSTERCOLLECTION_H__

#include <TNamed.h>
#include <vector>
#include "Cluster.h"

typedef std::map<Int_t,Cluster*> CL_CSET;

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// ClusterCollection                                                          //
//                                                                            //
// id[0]: unique id                                                           //
// id[1]: split cluster id (if exists)                                        // 
// id[2]: comprising cluster id (if exists)                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

class ClusterCollection: public TNamed {

public:

   ClusterCollection();
   ClusterCollection(Int_t uid);
   
   virtual ~ClusterCollection() {}
   
   virtual void clear() {}
   virtual void print(Int_t opt = 0) const;
  
   virtual void SetId(Int_t i, Int_t v); // !!! i > 0 !!! 
   
   //------------------------------------------------------------ 
 
   Int_t  GetId() const        { return Id_[0]; }
   Int_t  GetId(Int_t i) const { return Id_[i]; }
   
   Int_t  GetNClusters() const { return Clusters_.size(); }
 
   const CL_CSET& GetClusters() const { return Clusters_; }
   
   const Cluster* GetCluster(Int_t id) const;

   //------------------------------------------------------------ 
   
   virtual Bool_t FindCluster(const Cluster* cl) const ;
   virtual Bool_t FindCluster(Int_t id) const ;
   
   //------------------------------------------------------------ 
   
   virtual Bool_t AddCluster(Cluster* cl);
   
   virtual Bool_t RemoveCluster(Cluster* cl);
   virtual Bool_t RemoveCluster(Int_t id);
   virtual void   RemoveClusters();
  
protected:
  
    TArrayI Id_;
    CL_CSET Clusters_;  // <cluster ID, cluster>

    ClassDef(ClusterCollection,1)
};

#endif  /* __CLUSTERCOLLECTION_H__ */

