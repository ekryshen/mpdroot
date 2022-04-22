// $Id$
// Author: artur   2016/04/07

#ifndef __CLMANAGER_H__
#define __CLMANAGER_H__

#include <TObject.h>
#include <TString.h>
#include <TArrayD.h>

#include "ClusterCollection.h"

typedef std::map<Int_t, ClusterCollection *> CL_SSET;

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// ClManager                                                                  //
//                                                                            //
// <brief class description>                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

class ClManager {

public:
   ClManager();
   virtual ~ClManager();

   virtual void Print(Int_t opt = 0) const;
   virtual void DeleteAll(); // delete all elements, clusters and collections
   virtual void ClearCLPars() { fCLPars.Set(0); }

   //------------------------------------------------------------

   virtual Bool_t   InitCLPars();
   virtual Int_t    GetCLParsN();
   virtual Double_t GetCLPar(Int_t i);
   virtual void     SetCLPars(unsigned int ncopy, ...);

   //------------------------------------------------------------

   inline Int_t GetElementsN() const { return ElementsN_; }
   inline Int_t GetClustersN() const { return ClustersN_; }
   inline Int_t GetCollectionsN() const { return CollectionsN_; }

   virtual ClusterElement    *FindElement(int id) const;
   virtual Cluster           *FindCluster(int id) const;
   virtual ClusterCollection *FindCollection(int id) const;
   virtual ClusterCollection *FindCollection(TString name) const;

   //------------------------------------------------------------

   virtual ClusterElement    *CreateNewElement();
   virtual Cluster           *CreateNewCluster();
   virtual ClusterCollection *CreateNewCollection();

   //------------------------------------------------------------

   /* [1]- remove element from all clusters (without deleting)
      [1]- remove clusters ids from the element */
   virtual Bool_t RemoveElement(int id);

   /* [1]- remove cluster from all collections (without deleting)
      [1]- remove collections ids from the cluster

      [2]- remove all elements from the cluster (without deleting)
      [2]- remove cluster's id from the elements  */
   virtual Bool_t RemoveCluster(int id);

   /* [2]- remove all clusters from the collection (without deleting)
      [2]- remove collection's id from clusters */
   virtual Bool_t RemoveCollection(int id);

   virtual void DeleteElement(int id);    // remove & delete element anywhere
   virtual void DeleteCluster(int id);    // remove & delete cluster anywhere
   virtual void DeleteCollection(int id); // remove & delete collection anywhere

   virtual void DeleteCollections(); // remove & delete all collections anywhere

protected:
   CL_ESET Elements_;    // all-elements collection <element ID, element>
   CL_CSET Clusters_;    // all-clusters collection <cluster ID, cluster>
   CL_SSET Collections_; // all-cluster-sets collection <collection ID, collection>

   // private:

   Int_t ElementsN_;
   Int_t ClustersN_;
   Int_t CollectionsN_;

   TArrayD fCLPars;

   ClassDef(ClManager, 1);
};

#endif /* __CLMANAGER_H__ */
