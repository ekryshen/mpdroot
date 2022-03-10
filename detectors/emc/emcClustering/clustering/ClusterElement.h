// $Id$
// Author: artur   2016/04/07

#ifndef __CLUSTERELEMENT_H__
#define __CLUSTERELEMENT_H__

#include <TArrayI.h>
#include <TArrayD.h>
#include <TNamed.h>
#include <set>

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// ClusterElement                                                             //
//                                                                            //
// simplest data memebers interpretation:                                     //
//                                                                            //
// id[0]: unique id                                                           //
// id[1], id[2]: position it 2D-map                                           //
//                                                                            //
// value[0]: deposit                                                          //
// value[1]: deposit error                                                    //
// value[2]: background level                                                 //
//                                                                            //
// ClustersID: clusters to which the point belongs                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

typedef std::set<Int_t> CL_IDS;

class ClusterElement: public TNamed {

public:
  
   ClusterElement();
   ClusterElement(Int_t uid);
   ClusterElement(Int_t uid, Int_t id1, Int_t id2);
   
   virtual ~ClusterElement() {}
   
   virtual void clear(); // clear own data only
   virtual void print(Int_t opt = 0) const; 
   
   //------------------------------------------------------------
   
   virtual void SetId(Int_t i, Int_t id); // !!! i > 0 !!!
  
   virtual void SetValue(Double_t v)  { Value_[0] = v;  } 
   virtual void AddValue(Double_t v)  { Value_[0] += v; } 
   virtual void SetValue(Int_t i, Double_t v);
   virtual void AddValue(Int_t i, Double_t v);
   
   //------------------------------------------------------------
   
   Int_t    GetId()   const         { return Id_[0];    }
   Double_t GetValue()  const       { return Value_[0]; }  
   Int_t    GetId(Int_t i) const    { return Id_[i];    }
   Double_t GetValue(Int_t i) const { return Value_[i]; }
   
   const CL_IDS& GetClsId() const { return ClustersID_; }
   
   Int_t  GetNClusters()  const;
   
   Int_t* GetClustersID() const; 
   
   //------------------------------------------------------------
   
   virtual Bool_t FindCluster(Int_t id) const;
 
protected:
  
   TArrayI Id_;    // id
   TArrayD Value_; // own data
   
   CL_IDS  ClustersID_; // <cluster ID>
   
   /***** used in Cluster *****/
   
   Bool_t  AddCluster(Int_t id);
   Bool_t  RemoveCluster(Int_t id);
   
   friend class Cluster;

   ClassDef(ClusterElement,1)
};

#endif  /* __CLUSTERELEMENT_H__ */

