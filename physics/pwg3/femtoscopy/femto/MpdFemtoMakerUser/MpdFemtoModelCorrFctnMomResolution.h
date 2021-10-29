#ifndef ALIFEMTOMODELCORRFCTN3DKKGR_H
#define ALIFEMTOMODELCORRFCTN3DKKGR_H

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseCorrFctn.h"
#include "MpdFemtoModelManager.h"
#include "MpdFemtoModelHiddenInfo.h"

// ROOT
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH2D.h>
#include <TH3D.h>

class MpdFemtoModelCorrFctnMomResolution : public MpdFemtoBaseCorrFctn {
public:
    MpdFemtoModelCorrFctnMomResolution();
    MpdFemtoModelCorrFctnMomResolution(const char *title, Int_t aNbins, Double_t aQinvLo, Double_t aQinvHi);
    MpdFemtoModelCorrFctnMomResolution(const MpdFemtoModelCorrFctnMomResolution& aCorrFctn);
    virtual ~MpdFemtoModelCorrFctnMomResolution();

    MpdFemtoModelCorrFctnMomResolution& operator=(const MpdFemtoModelCorrFctnMomResolution& aCorrFctn);

    /// Method that allows front-loading model manager
    virtual void connectToManager(MpdFemtoModelManager *aManager);

    /// Make report
    virtual MpdFemtoString report();

    /// Add real pair
    virtual void addRealPair(MpdFemtoPair* aPair);
    /// Add mixed pair
    virtual void addMixedPair(MpdFemtoPair* aPair);

    /// Begin event
    virtual void eventBegin(const MpdFemtoEvent* aEvent);
    /// Event end
    virtual void eventEnd(const MpdFemtoEvent* aEvent);
    /// Finish
    virtual void finish();

    virtual void writeOutHistos();
    virtual TList* getOutputList();

    virtual MpdFemtoModelCorrFctnMomResolution* clone() const {
        return new MpdFemtoModelCorrFctnMomResolution(*this);
    }

    void SetFillkT(bool fillkT) {
        fFillkT = fillkT;
    }

    // void SetSpecificPairCut(AliFemtoPairCut* aCut);

    Double_t GetQinvTrue(MpdFemtoPair*);

    Double_t GetQoutTrue(MpdFemtoPair*);
    Double_t GetQsideTrue(MpdFemtoPair*);
    Double_t GetQlongTrue(MpdFemtoPair*);

    //Special MC analysis for K selected by PDG code -->
    void SetKaonPDG(Bool_t aSetKaonAna);

protected:
    MpdFemtoModelManager *fManager; // Link back to the manager to get the weights

    TH3D *fNumeratorTrue; // Numerator made with pairs from the same event
    TH3D *fNumeratorFake; // Numerator made with pairs from different events (mixed pairs)
    TH3D *fDenominator; // Denominator made with mixed pairs

    TH3D *fNumeratorTrueIdeal; // Numerator made with pairs (true qinv) from the same event
    TH3D *fNumeratorFakeIdeal; // Numerator made with pairs (true qinv) from different events (mixed pairs)
    TH3D *fDenominatorIdeal; // Denominator made with mixed pairs (true qinv)

    TH2D *fQgenQrec; // Qinv true (generated) vs. Qinv reconstructed


private:

    //Special MC analysis for K selected by PDG code -->
    Bool_t fKaonPDG;

    Bool_t fFillkT;
    const Int_t fNbbPairs = 21;
    TH1D *fkTdists[21]; // histograms with kT distributions for different BB pairs
    Double_t GetParentsKt(MpdFemtoPair *pair);
    Int_t GetPairNumber(MpdFemtoPair *pair); // returns pair code

    ClassDef(MpdFemtoModelCorrFctnMomResolution, 1);
};

#endif
