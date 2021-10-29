#include "MpdFemtoModelCorrFctnMomResolution.h"
//_______________________

MpdFemtoModelCorrFctnMomResolution::MpdFemtoModelCorrFctnMomResolution() :
fManager(nullptr),
fNumeratorTrue(nullptr),
fNumeratorFake(nullptr),
fDenominator(nullptr),
fNumeratorTrueIdeal(nullptr),
fNumeratorFakeIdeal(nullptr),
fDenominatorIdeal(nullptr),
fQgenQrec(nullptr),
fKaonPDG(kFALSE),
fFillkT(kFALSE) {
    // Default constructor
    fNumeratorTrue = new TH3D("ModelNumTrue", "ModelNumTrue", 50, 0.0, 0.5, 50, 0.0, 0.5, 50, 0.0, 0.5);
    fNumeratorFake = new TH3D("ModelNumFake", "ModelNumFake", 50, 0.0, 0.5, 50, 0.0, 0.5, 50, 0.0, 0.5);
    fDenominator = new TH3D("ModelDen", "ModelDen", 50, 0.0, 0.5, 50, 0.0, 0.5, 50, 0.0, 0.5);

    fNumeratorTrueIdeal = new TH3D("ModelNumTrueIdeal", "ModelNumTrueIdeal", 50, 0.0, 0.5, 50, 0.0, 0.5, 50, 0.0, 0.5);
    fNumeratorFakeIdeal = new TH3D("ModelNumFakeIdeal", "ModelNumFakeIdeal", 50, 0.0, 0.5, 50, 0.0, 0.5, 50, 0.0, 0.5);
    fDenominatorIdeal = new TH3D("ModelDenIdeal", "ModelDenIdeal", 50, 0.0, 0.5, 50, 0.0, 0.5, 50, 0.0, 0.5);

    fQgenQrec = new TH2D("QgenQrec", "QgenQrec", 50, 0.0, 0.5, 50, 0.0, 0.5);

    for (int i = 0; i < fNbbPairs; i++) {
        auto *name = Form("fkTdists[%i]", i);
        fkTdists[i] = new TH1D(name, name, 100, 0.0, 5.0);
        fkTdists[i]->Sumw2();
    }

    fNumeratorTrue->Sumw2();
    fNumeratorFake->Sumw2();
    fDenominator->Sumw2();

    fNumeratorTrueIdeal->Sumw2();
    fNumeratorFakeIdeal->Sumw2();
    fDenominatorIdeal->Sumw2();

    fQgenQrec->Sumw2();

}
//_______________________

MpdFemtoModelCorrFctnMomResolution::MpdFemtoModelCorrFctnMomResolution(const char *title,
        Int_t aNbins,
        Double_t aQinvLo,
        Double_t aQinvHi) :
fManager(nullptr),
fNumeratorTrue(nullptr),
fNumeratorFake(nullptr),
fDenominator(nullptr),
fNumeratorTrueIdeal(nullptr),
fNumeratorFakeIdeal(nullptr),
fDenominatorIdeal(nullptr),
fQgenQrec(nullptr),
fKaonPDG(kFALSE),
fFillkT(kFALSE) {
    // Normal constructor
    char *buf;
    buf = Form("NumTrue%s", title);
    fNumeratorTrue = new TH3D(buf, buf, aNbins, aQinvLo, aQinvHi, aNbins, aQinvLo, aQinvHi, aNbins, aQinvLo, aQinvHi);

    buf = Form("NumFake%s", title);
    fNumeratorFake = new TH3D(buf, buf, aNbins, aQinvLo, aQinvHi, aNbins, aQinvLo, aQinvHi, aNbins, aQinvLo, aQinvHi);

    buf = Form("Den%s", title);
    fDenominator = new TH3D(buf, buf, aNbins, aQinvLo, aQinvHi, aNbins, aQinvLo, aQinvHi, aNbins, aQinvLo, aQinvHi);

    buf = Form("NumTrueIdeal%s", title);
    fNumeratorTrueIdeal = new TH3D(buf, buf, aNbins, aQinvLo, aQinvHi, aNbins, aQinvLo, aQinvHi, aNbins, aQinvLo, aQinvHi);

    buf = Form("NumFakeIdeal%s", title);
    fNumeratorFakeIdeal = new TH3D(buf, buf, aNbins, aQinvLo, aQinvHi, aNbins, aQinvLo, aQinvHi, aNbins, aQinvLo, aQinvHi);

    buf = Form("DenIdeal%s", title);
    fDenominatorIdeal = new TH3D(buf, buf, aNbins, aQinvLo, aQinvHi, aNbins, aQinvLo, aQinvHi, aNbins, aQinvLo, aQinvHi);

    buf = Form("QgenQrec%s", title);
    fQgenQrec = new TH2D(buf, buf, aNbins, aQinvLo, aQinvHi, aNbins, aQinvLo, aQinvHi);
    //test
    //fQgenQrec = new TH2D(buf,buf,aNbins,aQinvLo,aQinvHi,aNbins,-0.05,0.05);

    for (int i = 0; i < fNbbPairs; i++) {
        const char *name = Form("fkTdists[%i]_%s", i, title);
        fkTdists[i] = new TH1D(name, name, 100, 0.0, 5.0);
        fkTdists[i]->Sumw2();
    }

    fNumeratorTrue->Sumw2();
    fNumeratorFake->Sumw2();
    fDenominator->Sumw2();

    fNumeratorTrueIdeal->Sumw2();
    fNumeratorFakeIdeal->Sumw2();
    fDenominatorIdeal->Sumw2();

    fQgenQrec->Sumw2();
}

//_______________________

MpdFemtoModelCorrFctnMomResolution::MpdFemtoModelCorrFctnMomResolution(const MpdFemtoModelCorrFctnMomResolution& aCorrFctn) :
fManager(aCorrFctn.fManager),
fNumeratorTrue(nullptr),
fNumeratorFake(nullptr),
fDenominator(nullptr),
fNumeratorTrueIdeal(nullptr),
fNumeratorFakeIdeal(nullptr),
fDenominatorIdeal(nullptr),
fQgenQrec(nullptr),
fKaonPDG(aCorrFctn.fKaonPDG),
fFillkT(aCorrFctn.fFillkT) {
    // Copy constructor
    fNumeratorTrue = new TH3D(*aCorrFctn.fNumeratorTrue);
    fNumeratorFake = new TH3D(*aCorrFctn.fNumeratorFake);
    fDenominator = new TH3D(*aCorrFctn.fDenominator);

    fNumeratorTrueIdeal = new TH3D(*aCorrFctn.fNumeratorTrueIdeal);
    fNumeratorFakeIdeal = new TH3D(*aCorrFctn.fNumeratorFakeIdeal);
    fDenominatorIdeal = new TH3D(*aCorrFctn.fDenominatorIdeal);

    fQgenQrec = new TH2D(*aCorrFctn.fQgenQrec);

    for (int i = 0; i < fNbbPairs; i++) {
        fkTdists[i] = new TH1D(*aCorrFctn.fkTdists[i]);
    }
}

//_______________________

MpdFemtoModelCorrFctnMomResolution::~MpdFemtoModelCorrFctnMomResolution() {
    // Destructor
    delete fNumeratorTrue;
    delete fNumeratorFake;
    delete fDenominator;

    delete fNumeratorTrueIdeal;
    delete fNumeratorFakeIdeal;
    delete fDenominatorIdeal;

    delete fQgenQrec;

    for (int i = 0; i < fNbbPairs; i++) {
        delete fkTdists[i];
    }

}
//_______________________

MpdFemtoModelCorrFctnMomResolution& MpdFemtoModelCorrFctnMomResolution::operator=(const MpdFemtoModelCorrFctnMomResolution& aCorrFctn) {
    // Assignment operator
    if (this == &aCorrFctn) {
        return *this;
    }

    MpdFemtoModelCorrFctnMomResolution::operator=(aCorrFctn);

    *fNumeratorTrue = *aCorrFctn.fNumeratorTrue;
    *fNumeratorFake = *aCorrFctn.fNumeratorFake;
    *fDenominator = *aCorrFctn.fDenominator;
    *fQgenQrec = *aCorrFctn.fQgenQrec;

    for (int i = 0; i < fNbbPairs; i++) {
        *fkTdists[i] = *aCorrFctn.fkTdists[i];
    }

    *fNumeratorTrueIdeal = *aCorrFctn.fNumeratorTrueIdeal;
    *fNumeratorFakeIdeal = *aCorrFctn.fNumeratorFakeIdeal;
    *fDenominatorIdeal = *aCorrFctn.fDenominatorIdeal;

    fManager = aCorrFctn.fManager;
    fKaonPDG = aCorrFctn.fKaonPDG;

    return *this;
}
//_______________________

void MpdFemtoModelCorrFctnMomResolution::connectToManager(MpdFemtoModelManager *aManager) {
    fManager = aManager;
}

//_______________________

MpdFemtoString MpdFemtoModelCorrFctnMomResolution::report() {
    // Construct the report
    TString report = "";

    return MpdFemtoString((const char *) report);
}

//_______________________

void MpdFemtoModelCorrFctnMomResolution::addRealPair(MpdFemtoPair* aPair) {
    if (mPairCut && !mPairCut->pass(aPair)) {
        return;
    }

    // if (!fKaonPDG) {
    // cout<<" AliFemtoModelCorrFcn add real pair "<<endl;
    Double_t weight = fManager->weight(aPair);
    //cout<<" wight "<< weight<<endl;
    //cout<<"Qinv"<<aPair->QInv()<<endl;

    fNumeratorTrue->Fill(aPair->qOutCMS(), aPair->qSideCMS(), aPair->qLongCMS(), weight); //qOut, qSide, qLong

    //Double_t tQinvTrue = GetQinvTrue(aPair);

    fNumeratorTrueIdeal->Fill(GetQoutTrue(aPair), GetQsideTrue(aPair), GetQlongTrue(aPair), weight); //qOut_ideal, qSide_ideal, qLong_ideal

    //cout<<"Qinv true"<<tQinvTrue<<endl;
    // }//Special MC analysis for K selected by PDG code -->
    //    else {
    //        Double_t weight = fManager->GetWeight(aPair);
    //        fNumeratorTrue->Fill(aPair->QOutCMS(), aPair->QSideCMS(), aPair->QLongCMS(), weight); //qOut, qSide, qLong
    //        //Double_t tQinvTrue = GetQinvTrue(aPair);
    //        fNumeratorTrueIdeal->Fill(GetQoutTrue(aPair), GetQsideTrue(aPair), GetQlongTrue(aPair), weight); //qOut_ideal, qSide_ideal, qLong_ideal
    //    }
}
//_______________________

void MpdFemtoModelCorrFctnMomResolution::addMixedPair(MpdFemtoPair* aPair) {
    if (mPairCut && !mPairCut->pass(aPair)) {
        return;
    }

    Double_t weight = fManager->weight(aPair);
    Double_t qinv = aPair->qInv();
    Double_t qinv_ideal = GetQinvTrue(aPair);

    //  if (!fKaonPDG) {
    fNumeratorFake->Fill(aPair->qOutCMS(), aPair->qSideCMS(), aPair->qLongCMS(), weight); //qOut, qSide, qLong
    fDenominator->Fill(aPair->qOutCMS(), aPair->qSideCMS(), aPair->qLongCMS(), 1.0); //qOut, qSide, qLong

    fNumeratorFakeIdeal->Fill(GetQoutTrue(aPair), GetQsideTrue(aPair), GetQlongTrue(aPair), weight); //qOut_ideal, qSide_ideal, qLong_ideal
    fDenominatorIdeal->Fill(GetQoutTrue(aPair), GetQsideTrue(aPair), GetQlongTrue(aPair), 1.0); //qOut_ideal, qSide_ideal, qLong_ideal

    if (fFillkT) {
        int pairNumber = GetPairNumber(aPair);

        if (pairNumber >= 0) {
            if (fkTdists[pairNumber]) {
                fkTdists[pairNumber]->Fill(GetParentsKt(aPair));
            }
        }
    }

    fQgenQrec->Fill(qinv_ideal, qinv);
    //  }//Special MC analysis for K selected by PDG code -->
    //    else {
    //        // AliFemtoTrack *inf1 = (AliFemtoTrack *) aPair->Track1()->Track();
    //        // AliFemtoTrack *inf2 = (AliFemtoTrack *) aPair->Track2()->Track();
    //        // Double_t pdg1 = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetPDGPid();
    //        // Double_t pdg2 = ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetPDGPid();
    //        // if((aPair->KT())<0.5)cout<<" Corr Func  pdg1 "<<pdg1<<" pdg2 "<<pdg2<<" qinv "<<aPair->QInv()<< " w "<<weight<<endl;
    //        fNumeratorFake->Fill(aPair->QOutCMS(), aPair->QSideCMS(), aPair->QLongCMS(), weight); //qOut, qSide, qLong
    //        fDenominator->Fill(aPair->QOutCMS(), aPair->QSideCMS(), aPair->QLongCMS(), 1.0); //qOut, qSide, qLong
    //        //if(qinv_ideal>0) {
    //        fNumeratorFakeIdeal->Fill(GetQoutTrue(aPair), GetQsideTrue(aPair), GetQlongTrue(aPair), weight); //qOut_ideal, qSide_ideal, qLong_ideal
    //        fDenominatorIdeal->Fill(GetQoutTrue(aPair), GetQsideTrue(aPair), GetQlongTrue(aPair), 1.0); //qOut_ideal, qSide_ideal, qLong_ideal
    //        fQgenQrec->Fill(qinv_ideal, qinv);
    //        //}
    //        //test
    //        //if(tQinvTrue>0)fQgenQrec->Fill(tQinvTrue,tQinvTrue-aPair->QInv());
    //    }
}

//void MpdFemtoModelCorrFctnMomResolution::SetSpecificPairCut(AliFemtoPairCut* aCut)
//{
//  fPairCut = aCut;
//}

//_______________________

Double_t MpdFemtoModelCorrFctnMomResolution::GetQinvTrue(MpdFemtoPair* aPair) {
    //if(!fKaonPDG) {

    MpdFemtoParticle* first = aPair->track1();
    MpdFemtoParticle* second = aPair->track2();

    if (!first || !second)
        return -1;

    MpdFemtoModelHiddenInfo* inf1 = (MpdFemtoModelHiddenInfo*) first->getHiddenInfo();
    MpdFemtoModelHiddenInfo* inf2 = (MpdFemtoModelHiddenInfo*) second->getHiddenInfo();

    if (!inf1 || !inf2)
        return -1;

    TVector3 mom1 = inf1->trueMomentum();
    Double_t mass1 = inf1->mass();
    Double_t ene1 = TMath::Sqrt(mom1.Mag2() + mass1 * mass1);

    TLorentzVector fm1(mom1.Px(), mom1.Py(), mom1.Pz(), ene1);

    TVector3 mom2 = inf2->trueMomentum();
    Double_t mass2 = inf2->mass();
    Double_t ene2 = TMath::Sqrt(mom2.Mag2() + mass2 * mass2);

    TLorentzVector fm2(mom2.Px(), mom2.Py(), mom2.Pz(), ene2);

    //std::cout<<" CFModel mass1 mass2 "<<am1<<" "<<am2<<std::endl;

    TLorentzVector tQinvTrueVec = fm1 - fm2;
    Double_t tQinvTrue = -1. * tQinvTrueVec.M();

    return tQinvTrue;
    //}
    //Special MC analysis for K selected by PDG code -->
    //  else {
    //    const AliFemtoTrack *inf1 = aPair->Track1()->Track(),
    //                        *inf2 = aPair->Track2()->Track();
    //
    //  AliFemtoLorentzVector fm1;
    //  AliFemtoThreeVector* temp = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetTrueMomentum();
    //  fm1.SetVect(*temp);
    //  Double_t am1 = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetMass();
    //  Double_t am2 = ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetMass();
    //
    //  am1=0.493677;
    //  am2=0.493677;
    //
    //  //Double_t pdg1 = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetPDGPid();
    //  //Double_t pdg2 = ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetPDGPid();
    //
    //  double ener = TMath::Sqrt(temp->Mag2()+am1*am1);
    //  fm1.SetE(ener);
    //
    //  AliFemtoLorentzVector fm2;
    //  AliFemtoThreeVector* temp2 =  ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetTrueMomentum();
    //  fm2.SetVect(*temp2);
    //  ener = TMath::Sqrt(temp2->Mag2()+am2*am2);
    //  fm2.SetE(ener);
    //
    //
    //  AliFemtoLorentzVector tQinvTrueVec = (fm1-fm2);
    //  Double_t tQinvTrue = -1.* tQinvTrueVec.m();
    //
    // // if(tQinvTrue<0 && am1!=0 && am2!=0)std::cout<<" CFModel Qinv mass1 mass2 "<<aPair->QInv()<<" Qinv_true "<<tQinvTrue<<" "<<am1<<" "<<am2<<" pdg1 "<<pdg1<<" pdg2 "<<pdg2<<std::endl;
    // // if(pdg1!=211 || pdg2!=211)std::cout<<" CFModel Qinv mass1 mass2 "<<aPair->QInv()<<" Qinv_true "<<tQinvTrue<<" "<<am1<<" "<<am2<<" pdg1 "<<pdg1<<" pdg2 "<<pdg2<<std::endl;
    //
    //  if(am1==0 || am2==0)tQinvTrue=-10000;
    //
    //  return tQinvTrue;
    //  }
}

//_______________________

Double_t MpdFemtoModelCorrFctnMomResolution::GetQoutTrue(MpdFemtoPair* aPair) {
    MpdFemtoParticle* first = aPair->track1();
    MpdFemtoParticle* second = aPair->track2();

    if (!first || !second)
        return -1;

    MpdFemtoModelHiddenInfo* inf1 = (MpdFemtoModelHiddenInfo*) first->getHiddenInfo();
    MpdFemtoModelHiddenInfo* inf2 = (MpdFemtoModelHiddenInfo*) second->getHiddenInfo();

    if (!inf1 || !inf2)
        return -1;

    TVector3 p1 = inf1->trueMomentum();
    TVector3 p2 = inf2->trueMomentum();

    Double_t dx = p1.Px() - p2.Px();
    Double_t dy = p1.Py() - p2.Py();

    Double_t px = p1.Px() + p2.Px();
    Double_t py = p1.Py() + p2.Py();

    Double_t qTkT = dx * px + dy * py;
    Double_t pT = TMath::Sqrt(px * px + py * py);
    // Double_t qOut;

    // if (!pT)
    //    qOut = 0;
    // else
    //    qOut = qTkT / pT;

    // if (!qOut)
    //    qOut = -10000;

    //return qOut;
    return qTkT / pT;
}

Double_t MpdFemtoModelCorrFctnMomResolution::GetQsideTrue(MpdFemtoPair* aPair) {
    MpdFemtoParticle* first = aPair->track1();
    MpdFemtoParticle* second = aPair->track2();

    if (!first || !second)
        return -1;

    MpdFemtoModelHiddenInfo* inf1 = (MpdFemtoModelHiddenInfo*) first->getHiddenInfo();
    MpdFemtoModelHiddenInfo* inf2 = (MpdFemtoModelHiddenInfo*) second->getHiddenInfo();

    if (!inf1 || !inf2)
        return -1;

    TVector3 p1 = inf1->trueMomentum();
    TVector3 p2 = inf2->trueMomentum();

    Double_t x1 = p1.Px();
    Double_t x2 = p2.Px();
    Double_t y1 = p1.Py();
    Double_t y2 = p2.Py();

    Double_t px = p1.Px() + p2.Px();
    Double_t py = p1.Py() + p2.Py();

    Double_t qTkT = x2 * y1 - x1 * y2;
    Double_t pT = sqrt(px * px + py * py);
    //   Double_t qSide;

    //    if (!pT)
    //        qSide = 0;
    //    else
    //        qSide = (2.0 * qTkT) / pT;
    //
    //    if (!qSide)
    //        qSide = -10000;
    //
    //    return qSide;
    return (2. * qTkT) / pT;
}

Double_t MpdFemtoModelCorrFctnMomResolution::GetQlongTrue(MpdFemtoPair* aPair) {
    MpdFemtoParticle* first = aPair->track1();
    MpdFemtoParticle* second = aPair->track2();

    if (!first || !second)
        return -1;

    MpdFemtoModelHiddenInfo* inf1 = (MpdFemtoModelHiddenInfo*) first->getHiddenInfo();
    MpdFemtoModelHiddenInfo* inf2 = (MpdFemtoModelHiddenInfo*) second->getHiddenInfo();

    if (!inf1 || !inf2)
        return -1;

    TVector3 mom1 = inf1->trueMomentum();
    Double_t mass1 = inf1->mass();
    Double_t ene1 = TMath::Sqrt(mom1.Mag2() + mass1 * mass1);

    TLorentzVector fm1(mom1.Px(), mom1.Py(), mom1.Pz(), ene1);

    TVector3 mom2 = inf2->trueMomentum();
    Double_t mass2 = inf2->mass();
    Double_t ene2 = TMath::Sqrt(mom2.Mag2() + mass2 * mass2);

    TLorentzVector fm2(mom2.Px(), mom2.Py(), mom2.Pz(), ene2);

    //std::cout<<" CFModel mass1 mass2 "<<am1<<" "<<am2<<std::endl;

    Double_t dz = fm1.Z() - fm2.Z();
    Double_t zz = fm1.Z() + fm2.Z();

    Double_t dt = fm1.T() - fm2.T();
    Double_t tt = fm1.T() + fm2.T();

    Double_t beta = zz / tt;
    Double_t gamma = 1.0 / TMath::Sqrt((1. - beta)*(1. + beta));

    Double_t qLong = gamma * (dz - beta * dt);

    //    if (!qLong)
    //        qLong = -10000;

    return qLong;
}

//_______________________

void MpdFemtoModelCorrFctnMomResolution::eventBegin(const MpdFemtoEvent* /* aEvent */) {
    /* Do nothing */
}
//_______________________

void MpdFemtoModelCorrFctnMomResolution::eventEnd(const MpdFemtoEvent* /* aEvent */) {
    /* Do nothing */
}
//_______________________

void MpdFemtoModelCorrFctnMomResolution::finish() {
    /* Do nothing */
}
//_______________________

void MpdFemtoModelCorrFctnMomResolution::writeOutHistos() {
    // Write out data histos

    fQgenQrec->Write();

    fNumeratorTrue->Write();
    fNumeratorFake->Write();
    fDenominator->Write();

    if (fFillkT) {
        for (int i = 0; i < fNbbPairs; i++) {
            fkTdists[i]->Write();
        }
    }
    fNumeratorTrueIdeal->Write();
    fNumeratorFakeIdeal->Write();
    fDenominatorIdeal->Write();


}
//_________________________

TList* MpdFemtoModelCorrFctnMomResolution::getOutputList() {
    // Prepare the list of objects to be written to the output
    TList *tOutputList = new TList();

    tOutputList->Add(fNumeratorTrue);
    tOutputList->Add(fNumeratorFake);
    tOutputList->Add(fDenominator);

    tOutputList->Add(fNumeratorTrueIdeal);
    tOutputList->Add(fNumeratorFakeIdeal);
    tOutputList->Add(fDenominatorIdeal);
    tOutputList->Add(fQgenQrec);

    if (fFillkT) {
        for (int i = 0; i < fNbbPairs; i++) {
            tOutputList->Add(fkTdists[i]);
        }
    }
    return tOutputList;
}

void MpdFemtoModelCorrFctnMomResolution::SetKaonPDG(Bool_t aSetKaonAna) {
    fKaonPDG = aSetKaonAna;
}

Double_t MpdFemtoModelCorrFctnMomResolution::GetParentsKt(MpdFemtoPair *pair) {
//    MpdFemtoParticle* first = new MpdFemtoParticle(*(pair->track1()));
//    MpdFemtoParticle* second = new MpdFemtoParticle(*(pair->track2()));
//
//    if (!first) {
//        if (second)
//            delete second;
//        return -1;
//    }
//    
//    if (!second) {
//        if (first)
//            delete first;
//        return -1;
//    }
//
//    MpdFemtoModelHiddenInfo* info1 = first->getHiddenInfo();
//    MpdFemtoModelHiddenInfo* info2 = second->getHiddenInfo();
//
//    if (!info1 || !info2) {
//        if (first) delete first;
//        if (second) delete second;
//        return -1;
//    }
//    AliFemtoThreeVector* p1 = info1->GetMotherMomentum();
//    AliFemtoThreeVector* p2 = info2->GetMotherMomentum();
//
//    if (!p1 || !p2) {
//        if (first) delete first;
//        if (second) delete second;
//        return -1;
//    }
//    double px = p1->x() + p2->x();
//    double py = p1->y() + p2->y();
//    double pT = sqrt(px * px + py * py);
//
//    delete first;
//    delete second;
//
//    return pT / 2.;
    return 0.;
}
//
Int_t MpdFemtoModelCorrFctnMomResolution::GetPairNumber(MpdFemtoPair *pair) {
//    const AliFemtoModelHiddenInfo
//            *info1 = (AliFemtoModelHiddenInfo*) pair->Track1()->GetHiddenInfo(),
//            *info2 = (AliFemtoModelHiddenInfo*) pair->Track2()->GetHiddenInfo();
//
//    if (!info1 || !info2) {
//        return -1;
//    }
//
//    int pdg1 = TMath::Abs(info1->GetMotherPdgCode());
//    int pdg2 = TMath::Abs(info2->GetMotherPdgCode());
//
//    if (pdg2 < pdg1) {
//        std::swap(pdg1, pdg2);
//    }
//
//    if (pdg1 == 2212 && pdg2 == 2212) return 0; // pp
//    if (pdg1 == 2212 && pdg2 == 3122) return 1; // pΛ
//    if (pdg1 == 3122 && pdg2 == 3122) return 2; // ΛΛ
//    if (pdg1 == 2212 && pdg2 == 3222) return 3; // pΣ+
//    if (pdg1 == 3122 && pdg2 == 3222) return 4; // ΛΣ+
//    if (pdg1 == 3222 && pdg2 == 3222) return 5; // Σ+Σ+
//    if (pdg1 == 2212 && pdg2 == 3312) return 6; // pΞ-
//    if (pdg1 == 2212 && pdg2 == 3322) return 7; // pΞ0
//    if (pdg1 == 3122 && pdg2 == 3312) return 8; // ΛΞ-
//    if (pdg1 == 3122 && pdg2 == 3322) return 9; // ΛΞ0
//    if (pdg1 == 3222 && pdg2 == 3322) return 10; // Σ+Ξ0
//    if (pdg1 == 3222 && pdg2 == 3312) return 11; // Σ+Ξ-
//    if (pdg1 == 3212 && pdg2 == 3122) return 12; // Σ0Λ
//    if (pdg1 == 2212 && pdg2 == 3212) return 13; // pΣ0
//    if (pdg1 == 3212 && pdg2 == 3222) return 14; // Σ0Σ+
//    if (pdg1 == 3322 && pdg2 == 3322) return 15; // Ξ0Ξ0
//    if (pdg1 == 3312 && pdg2 == 3322) return 16; // Ξ-Ξ0
//    if (pdg1 == 3312 && pdg2 == 3312) return 17; // Ξ-Ξ-
//    if (pdg1 == 3212 && pdg2 == 3322) return 18; // Σ0Ξ0
//    if (pdg1 == 3212 && pdg2 == 3312) return 19; // Σ0Ξ-
//    if (pdg1 == 3212 && pdg2 == 3212) return 20; // Σ0Σ0
//
//    return -1;
    return 0;
}
