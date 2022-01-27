/*************************************************************************************
 *
 *         Class MpdZdcDigiScheme
 *
 *  Author:   Elena Litvinenko
 *  e-mail:   litvin@nf.jinr.ru
 *  Version:  16-May-2008
 *
 ************************************************************************************/

#ifndef MPDZDCDIGISCHEME_H
#define MPDZDCDIGISCHEME_H

#include "FairGeoNode.h"
#include "MpdZdcGeoPar.h"
#include "MpdZdcDigiPar.h"
#include "MpdZdcPsdGeoPar.h"
#include "MpdZdcPsdDigiPar.h"

#include "TObject.h"
#include <map>
#include <vector>

using std::pair;
using std::vector;

typedef std::vector<Int_t>    MpdZdcVolId_t;   // now - {MotherMotherCopyNo, MotherCopyNo, VolumeId, CopyNo}
typedef std::vector<Int_t>    MpdZdcDigiId_t;  // now - {DetectorID, ModuleID, ChannelID}
typedef std::vector<Double_t> MpdZdcVolInfo_t; // now - Center X,Y,Z; and Dx,Dy,Dz from center to the corner

class MpdZdcDigiScheme : public TObject {
public:
   MpdZdcDigiScheme();
   virtual ~MpdZdcDigiScheme();

   Bool_t Init(MpdZdcGeoPar *geoPar, MpdZdcDigiPar *digiPar, Bool_t pAddPsd = kFALSE, Int_t pVerbose = 0);
   Bool_t InitPsd(MpdZdcPsdGeoPar *geoPar, MpdZdcPsdDigiPar *digiPar, Int_t pVerbose = 0);

   static MpdZdcDigiScheme *Instance();

   Bool_t AddNodes(TObjArray *sensNodes, Int_t pGlobalDetectorNumber, Bool_t pAddPsd = kFALSE, Int_t pVerbose = 0);
   MpdZdcVolInfo_t *CreateVolInfoElement(FairGeoNode *nod, Int_t pVerbose);
   MpdZdcVolId_t *  CreateVolElement(FairGeoNode *nod, Int_t nodeNumber, MpdZdcDigiId_t *right,
                                     Int_t pGlobalDetectorNumber, Int_t pVerbose);
   Bool_t           CreateVolCopyElements(MpdZdcVolId_t *left, MpdZdcDigiId_t *right);
   Bool_t           CreateVolInfoCopyElements(MpdZdcDigiId_t *right, MpdZdcVolInfo_t *volInfo);

   Bool_t         IsVolumeExist(MpdZdcVolId_t *pVolId);
   MpdZdcDigiId_t GetDigiId(MpdZdcVolId_t *pVolId);
   Int_t          GetDetectorID(MpdZdcVolId_t *pVolId);
   Int_t          GetChannelID(MpdZdcVolId_t *pVolId);
   Bool_t         GetDetIdModIdChanId(Int_t pMcVolumeNumber, Int_t pMcCopyNumber, Int_t pMotherCopyNumber,
                                      Int_t pMotherMotherCopyNumber, Int_t &pDetId, Int_t &pChanId, Int_t &pModId);

   Bool_t           GetVolCenterXYZ(MpdZdcDigiId_t *pDigiId, Double_t &x, Double_t &y, Double_t &z);
   Bool_t           GetVolDxDyDz(MpdZdcDigiId_t *pDigiId, Double_t &Dx, Double_t &Dy, Double_t &Dz);
   MpdZdcVolInfo_t *GetVolInfo(MpdZdcVolId_t *pVolId);

   void Print();
   void PrintVolume(Int_t volID, Int_t copyNo, Int_t copyNoMother, Int_t copyNoMotherMother = 1);

   Bool_t CalcDimensions(Int_t pGlobalDetectorNumber, Int_t &nx, Int_t &ny, Int_t &nz);
   void   GetZdcDimensions(Int_t &nx, Int_t &ny, Int_t &nz);
   void   GetZdcPsdDimensions(Int_t &nx_psd, Int_t &ny_psd, Int_t &nz_psd);

   MpdZdcDigiId_t GetDigiIdFromCoords(Double_t x, Double_t y, Double_t z);
   MpdZdcDigiId_t GetDigiIdFromVolumeData(Int_t pMcVolumeNumber, Int_t pMcCopyNumber, Int_t pMotherCopyNumber,
                                          Int_t pMotherMotherCopyNumber);
   void           SplitDigiID(MpdZdcDigiId_t digiID, Int_t &detID, Int_t &modID, Int_t &chanID);

   inline MpdZdcDigiPar *   GetZdcDigiPar() { return fZdcDigiPar; };
   inline MpdZdcPsdDigiPar *GetZdcPsdDigiPar() { return fZdcPsdDigiPar; };

protected:
   static MpdZdcDigiScheme *fInstance;    // Instance of singleton object
   static Int_t             fRefcount;    // Counter of references on this
   static Bool_t            fInitialized; // Defines whether was initialized
   static Bool_t            fYesPsd;      //

private:
   std::map<MpdZdcVolId_t, MpdZdcDigiId_t> fVolToDigiIdMap; //! correspondence for all active volumes (Zdc+ZdcPsd)
   std::map<MpdZdcDigiId_t, MpdZdcVolInfo_t *>
      fDigiToVolInfoMap; //! correspondence for active volumes (must exist for ZdcPsd, can also exist for Zdc)

   Int_t Nx;
   Int_t Ny;
   Int_t Nz;
   Int_t Nx_psd;
   Int_t Ny_psd;
   Int_t Nz_psd;

   MpdZdcDigiPar *   fZdcDigiPar;    //!
   MpdZdcPsdDigiPar *fZdcPsdDigiPar; //!
   TObjArray *       fPasNodes;      //!

   ClassDef(MpdZdcDigiScheme, 1);
};
#endif // MPDZDCDIGISCHEME_H
