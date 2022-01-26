// $Id$
// Author: artur   2016/04/15

#ifndef __EMCCLMANAGER_H__
#define __EMCCLMANAGER_H__

#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include "FairTask.h"
#include "ClManager.h"
#include "EmcCLinfo.h"
#include "EmcClusterCollection.h"
#include "MpdEmcHitA.h"
#include "MpdEmcGeoParWrapper.h"

class TClonesArray;
class TObjArray;
class EmcClusterInfo;

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// EmcClManager                                                               //
//                                                                            //
// <brief class description>                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

class EmcClManager: public FairTask,  public ClManager {

public:

    EmcClManager();
    virtual ~EmcClManager();
    
    virtual InitStatus Init();   
    virtual void Exec(Option_t* opt);
    virtual void Finish();
    
    virtual EmcClusterElement*    CreateNewElement();
    virtual EmcCluster*           CreateNewCluster();
    virtual EmcClusterCollection* CreateNewCollection();
    
    void   ResetCounter() { fNum = 0; }
    
    virtual void Clear();
    
    Bool_t SearchHitsForErrors(TClonesArray* digits = 0);
    TH2D*  CreateHitMapHisto(TString name, TString title) const;
    
    virtual void Fill(TClonesArray* digits = 0);
    
    virtual void FillClusterInfo(Int_t id, TClonesArray* clinfo, 
                                 Int_t clflag, TH2D* htot = 0, TH2D* hcl = 0);  
    virtual void FillClusterInfo(TString coll_name, TClonesArray* clinfo, 
                                 Int_t clflag, TH2D* htot = 0, TH2D* hcl = 0);
    virtual void FillClusterInfo(EmcClusterCollection* clcoll, TClonesArray* clinfo, 
                                 Int_t clflag, TH2D* htot = 0, TH2D* hcl = 0);
    
    virtual void FillClusterInfoHisto(TClonesArray* clinfo, TH2D* hcl);
    
    inline TH2D* GetHitsMap()          { return fHitsMap; }
    inline TH2D* GetHitsMapSelected()  { return fHitsMapSelected; }
    inline TH2D* GetHitsMapClustered() { return fHitsMapClustered; }
    
    inline TH2D* GetMCClusters()       { return fMCClusters; }
    inline TH2D* GetRecoClusters()     { return fRCClusters; }
    
    inline TClonesArray*  GetMCClustersInfo()   { return fMCClustersInfo; }
    inline TClonesArray*  GetRecoClustersInfo() { return fRCClustersInfo; }
       
    void SaveCLinfo(TString filename, TString path = "");
    
    virtual void Print(Int_t opt = 0) const;
    
protected:
  
    EmcCLinfo*           fCLinfo;
    MpdEmcGeoParWrapper* fGeoPar;
    TClonesArray*        fHitsArray;
    
    TClonesArray*        fMCClustersInfo;
    TClonesArray*        fRCClustersInfo;
    
    TH2D*  fHitsMap;
    TH2D*  fHitsMapSelected;
    TH2D*  fHitsMapClustered;
    
    TH2D*  fMCClusters;
    TH2D*  fRCClusters;
    
    Int_t  fNum;
    
    Int_t  fCLMotherPdg;
    
    virtual void MakeMCClusters();
    virtual void RecoClusters_1(Float_t thr, Int_t ncut);
    
    void SaveHistos();
    
    ClassDef(EmcClManager,1)
};

#endif  /* __EMCCLMANAGER_H__ */

