// $Id$
// Author: artur   2016/04/07

#ifndef __CLUSTER_H__
#define __CLUSTER_H__

#include <TNamed.h>
#include <map>
#include <set>
#include "ClusterElement.h"

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Cluster                                                                    //
//                                                                            //
// simplest data memebers interpretation:                                     //
//                                                                            //
// id[0]: unique id                                                           //
// id[1], id[2]: mother particle & track number                               //
//                                                                            //
// value[0]: deposit                                                          //
// value[1], value[2]: exact position it 2D-map                               //
//                                                                            //
// CollectionsID: collections to which the cluster belongs                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

typedef std::map<Int_t,ClusterElement*> CL_ESET;

class Cluster: public TNamed {
  
public:
  
   Cluster();
   Cluster(Int_t uid);
   
   virtual ~Cluster() {}
   
   virtual void clear(); // clear own data only 
   virtual void print(Int_t opt = 0) const;
     
   //------------------------------------------------------------ 
   
   virtual void SetId(Int_t i, Int_t id); // !!! i > 0 !!!
   
   virtual void SetValue(Double_t v) { Value_[0] = v;  } 
   virtual void AddValue(Double_t v) { Value_[0] += v; } 
   virtual void SetValue(Int_t i, Double_t v);
   virtual void AddValue(Int_t i, Double_t v);
   
   //------------------------------------------------------------
   
   Int_t    GetId()  const            { return Id_[0];    }
   Double_t GetValue() const          { return Value_[0]; }  
   Int_t    GetId(Int_t i) const      { return Id_[i];    }
   Double_t GetValue(Int_t i) const   { return Value_[i]; }
      
   const CL_ESET& GetElements() const { return Elements_; }
   const CL_IDS&  GetCollsID()  const { return CollectionsID_; }
   
   Int_t  GetNElements()    const { return Elements_.size(); }
   Int_t  GetNCollections() const { return CollectionsID_.size(); }
   
   Int_t* GetElementsID()    const; 
   Int_t* GetCollectionsID() const; 
   
   const ClusterElement* GetElement(Int_t id) const;
   
   //------------------------------------------------------------ 
  
   virtual Bool_t   FindElement(const ClusterElement* elem) const;
   virtual Bool_t   FindElement(Int_t id) const;
   virtual Bool_t   FindCollection(Int_t id) const;
   
   virtual Bool_t   AddElement(ClusterElement* elem);
   
   virtual Bool_t   RemoveElement(ClusterElement* elem);
   virtual Bool_t   RemoveElement(Int_t id);
   virtual void     RemoveElements();
    
protected:
   
    TArrayI  Id_;    // id 
    TArrayD  Value_; // own data
    
    CL_ESET  Elements_;       // <element ID, element>
    CL_IDS   CollectionsID_;  // <collection ID>
    
    /***** used in ClusterCollection *****/
    
    virtual Bool_t AddCollection(Int_t id);
    virtual Bool_t RemoveCollection(Int_t id);
    
    friend class ClusterCollection;
       
    ClassDef(Cluster,1)
};

#endif  /* __CLUSTER_H__ */

