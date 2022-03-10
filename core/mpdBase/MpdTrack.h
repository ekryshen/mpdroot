// Author: Emelyanov D.
// Update: Oleg Rogachevsky 2009-09-17 17:42:25+0400
// Copyright: 2009 (C) Oleg Rogachevsky

#ifndef ROOT_MpdTrack
#define ROOT_MpdTrack
#ifndef ROOT_TObject
#include "TObject.h"
#endif
#include "MpdHelix.h"

class MpdTrack : public TObject {

 private:

    Int_t     fID; // Track number from reconstruction  
    Int_t     fNofHits; // Number of track hits
    Int_t     fNofHitsPossTpc; //   
    Int_t     fNofHitsFitTpc; // 
    
    Float_t   fPidProbElectron;  // [combined] probability  
    Float_t   fPidProbPion;  // [combined] probability  
    Float_t   fPidProbKaon;  // [combined] probability
    Float_t   fPidProbProton;  // [combined] probability
    
    Float_t   fPidTPCProbElectron;  // [tpc] probability  
    Float_t   fPidTPCProbPion;  // [tpc] probability  
    Float_t   fPidTPCProbKaon;  // [tpc] probability
    Float_t   fPidTPCProbProton;  // [tpc] probability
    
    Float_t   fPidTOFProbElectron;  // [tof] probability  
    Float_t   fPidTOFProbPion;  // [tof] probability  
    Float_t   fPidTOFProbKaon;  // [tof] probability
    Float_t   fPidTOFProbProton;  // [tof] probability
       
    Float_t   fTofBeta;   // [tof] beta
    Float_t   fTofMass2;  // [tof] identified m^2, can be >0, <0 and Nan!
    Float_t   fdEdXTPC;  // [tof]
    Int_t     fTofHitIndex;  // [tof] 
    Int_t     fTofFlag;  // [tof] it is set if tof identification exists
    Float_t   fNSigmaElectron;  //   
    Float_t   fNSigmaPion;  //   
    Float_t   fNSigmaKaon;  //   
    Float_t   fNSigmaProton;  //   
    Float_t   fChi2;  // Track Chi2
    Float_t   fPt;  // Signed Transverce momentum  
    Float_t   fTheta;  // Theta angle (from beam line)  
    Float_t   fPhi;  // Phi angle  
    Float_t   fPtError;  // Pt error   
    Float_t   fThetaError;  // Theta error  
    Float_t   fPhiError;  // Phi error  
    Float_t   fDCAX;  //   
    Float_t   fDCAY;  //   
    Float_t   fDCAZ;  //   
    Float_t   fDCAGlobalX;  //   
    Float_t   fDCAGlobalY;  //   
    Float_t   fDCAGlobalZ;  //   
    Float_t   fFirstPointX;  // point closest to beam line
    Float_t   fFirstPointY;  // point closest to beam line
    Float_t   fFirstPointZ;  // point closest to beam line
    Float_t   fLastPointX;  //   
    Float_t   fLastPointY;  //   
    Float_t   fLastPointZ;  //   
    Short_t   fHelixQ;  //   
    Bool_t    fEdgeCut;  // kTRUE if number of hits closer to boundaries than 1.5 cm divided by nHits is larger than 50% (else: kFALSE)
    ULong64_t fHitMap; // n-th bit is 1 if there is hit in layer, 0 otherwise
    ULong64_t fSharedHitMap;// n-th bit if bit if hit in given layer shared with any other track

 public:

  MpdTrack(); // Default constructor
  virtual ~MpdTrack() {} // Destructor

    void SetID( Int_t n ) {fID=n;}
    void SetNofHits( Int_t n ) {fNofHits=n;}
    void SetNofHitsPossTpc( Int_t n ) {fNofHitsPossTpc=n;}
    void SetNofHitsFitTpc( Int_t n ) {fNofHitsFitTpc=n;}
    void SetPidProbElectron( Float_t n ) {fPidProbElectron=n;} 
    void SetPidProbPion( Float_t n ) {fPidProbPion=n;} 
    void SetPidProbKaon( Float_t n ) {fPidProbKaon=n;} 
    void SetPidProbProton( Float_t n ) {fPidProbProton=n;} 
    void SetTOFpidProb( Float_t n1, Float_t n2, Float_t n3, Float_t n4, Int_t flag ) {
    	fPidTOFProbElectron=n1; fPidTOFProbPion=n2; fPidTOFProbKaon=n3; fPidTOFProbProton=n4; SetTofFlag(flag);} 
    void SetTPCpidProb( Float_t n1, Float_t n2, Float_t n3, Float_t n4,  Int_t flag) {
    	fPidTPCProbElectron=n1; fPidTPCProbPion=n2; fPidTPCProbKaon=n3; fPidTPCProbProton=n4; SetTofFlag(flag);}     
    void SetCombPidProb( Float_t n1, Float_t n2, Float_t n3, Float_t n4) {
    	fPidProbElectron=n1; fPidProbPion=n2; fPidProbKaon=n3; fPidProbProton=n4;} 
    void SetTofBeta(Float_t n) {fTofBeta=n;}  
    void SetTofMass2(Float_t n) {fTofMass2=n;}  
    void SetdEdXTPC(Float_t n) {fdEdXTPC=n;}  
    void SetTofHitIndex(Int_t n) {fTofHitIndex=n;}  
    void SetTofFlag(Int_t n) {fTofFlag = fTofFlag | n;} 
    void SetNSigmaElectron( Float_t n ) {fNSigmaElectron=n;}
    void SetNSigmaPion( Float_t n ) {fNSigmaPion=n;}
    void SetNSigmaKaon( Float_t n ) {fNSigmaKaon=n;}
    void SetNSigmaProton( Float_t n ) {fNSigmaProton=n;}
    void SetChi2( Float_t n ) {fChi2=n;}
    void SetPt( Float_t n ) {fPt=n;}
    void SetTheta( Float_t n ) {fTheta=n;}
    void SetPhi( Float_t n ) {fPhi=n;}
    void SetPtError( Float_t n ) {fPtError=n;}
    void SetThetaError( Float_t n ) {fThetaError=n;}
    void SetPhiError( Float_t n ) {fPhiError=n;}
    void SetDCAX( Float_t n ) {fDCAX=n;}
    void SetDCAY( Float_t n ) {fDCAY=n;}
    void SetDCAZ( Float_t n ) {fDCAZ=n;}
    void SetDCAGlobalX( Float_t n ) {fDCAGlobalX=n;}
    void SetDCAGlobalY( Float_t n ) {fDCAGlobalY=n;}
    void SetDCAGlobalZ( Float_t n ) {fDCAGlobalZ=n;}
    void SetFirstPointX( Float_t n ) {fFirstPointX=n;}
    void SetFirstPointY( Float_t n ) {fFirstPointY=n;}
    void SetFirstPointZ( Float_t n ) {fFirstPointZ=n;}
    void SetLastPointX( Float_t n ) {fLastPointX=n;}
    void SetLastPointY( Float_t n ) {fLastPointY=n;}
    void SetLastPointZ( Float_t n ) {fLastPointZ=n;}
    void SetHelixQ( Short_t n ) {fHelixQ=n;}
    void SetEdgeCut( Bool_t n ) {fEdgeCut=n;}
    void SetLayerHitMap(ULong64_t map){fHitMap = map;};
    void SetSharedHitMap(ULong64_t map){fSharedHitMap = map;}

    Int_t     GetID()const {return fID;}
    Int_t     GetNofHits()const {return fNofHits;}
    Int_t     GetNofHitsPossTpc()const {return fNofHitsPossTpc;}
    Int_t     GetNofHitsFitTpc()const {return fNofHitsFitTpc;}
    
    Float_t   GetTPCPidProbElectron()const {return fPidTPCProbElectron;}
    Float_t   GetTPCPidProbPion()const {return fPidTPCProbPion;}
    Float_t   GetTPCPidProbKaon()const {return fPidTPCProbKaon;}
    Float_t   GetTPCPidProbProton()const {return fPidTPCProbProton;}
    
    Float_t   GetTOFPidProbElectron()const {return fPidTOFProbElectron;}
    Float_t   GetTOFPidProbPion()const {return fPidTOFProbPion;}
    Float_t   GetTOFPidProbKaon()const {return fPidTOFProbKaon;}
    Float_t   GetTOFPidProbProton()const {return fPidTOFProbProton;}
    
    Float_t   GetPidProbElectron()const {return fPidProbElectron;}
    Float_t   GetPidProbPion()const {return fPidProbPion;}
    Float_t   GetPidProbKaon()const {return fPidProbKaon;}
    Float_t   GetPidProbProton()const {return fPidProbProton;}
      
    Float_t   GetTofBeta() const{return fTofBeta;}
    Float_t   GetTofMass2() const{return fTofMass2;}
    Float_t   GetdEdXTPC() const{return fdEdXTPC;}
    Int_t     GetTofHitIndex() const{return fTofHitIndex;};
    Int_t     GetTofFlag() const{return fTofFlag;}
    Float_t   GetNSigmaElectron() const{return fNSigmaElectron;}
    Float_t   GetNSigmaPion() const{return fNSigmaPion;}
    Float_t   GetNSigmaKaon() const{return fNSigmaKaon;}
    Float_t   GetNSigmaProton() const{return fNSigmaProton;}
    Float_t   GetChi2() {return fChi2;}
    Float_t   GetPt() const{return fPt;}
    Float_t   GetTheta() const{return fTheta;}
    Float_t   GetPhi() const{return fPhi;}
    Float_t   GetPtError() const {return fPtError;}
    Float_t   GetThetaError()const {return fThetaError;}
    Float_t   GetPhiError() const{return fPhiError;}
    Float_t   GetPx() const;
    Float_t   GetPy() const;
    Float_t   GetPz() const;
    Float_t   GetEta() const;
    Float_t   GetDCAX() const{return fDCAX;}
    Float_t   GetDCAY() const {return fDCAY;}
    Float_t   GetDCAZ() const {return fDCAZ;}
    Float_t   GetDCAGlobalX() const{return fDCAGlobalX;}
    Float_t   GetDCAGlobalY() const {return fDCAGlobalY;}
    Float_t   GetDCAGlobalZ() const{return fDCAGlobalZ;}
    Float_t   GetFirstPointX() const {return fFirstPointX;}
    Float_t   GetFirstPointY() const {return fFirstPointY;}
    Float_t   GetFirstPointZ() const {return fFirstPointZ;}
    Float_t   GetLastPointX() const {return fLastPointX;}
    Float_t   GetLastPointY() const{return fLastPointY;}
    Float_t   GetLastPointZ() const{return fLastPointZ;}
    Float_t   GetCharge() const{return -TMath::Sign(1.0,GetPt());}
    Short_t   GetHelixQ() const{return fHelixQ;}
    Bool_t    GetEdgeCut() {return fEdgeCut;}
    ULong64_t GetLayerHitMap()const{return fHitMap;};
    ULong64_t GetSharedHitMap()const{return fSharedHitMap;};
    Int_t GetNSharedTpcHits()const;
    MpdHelix GetHelix()const;
  ClassDef(MpdTrack,4)
};

#endif
