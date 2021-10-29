// =======================================================================================================
// This code is based on the Identity module created by
// A. Rustamov and M. Arslandok, NIM A 946 (2019) 162622, https://arxiv.org/abs/1807.06370,
// and also a script provided by A.Rustamov.
// This version of the code is prepared by Igor Altsybeev.
// =======================================================================================================

#include "src/TIdentity2D.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVectorF.h"
#include "TSystem.h"
#include "TObject.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TMath.h"
#include <iomanip>
#include "iostream"
#include "string"

using namespace std;
using std::cout;
using std::setw;

const int nSubs = 10;

// =======================================================================================================
// Helper Functions
double calc_bcorr( double cov, double f, double b, double f2, double b2 )
{
    return (cov-f*b)/sqrt( (f2-f*f) * (b2-b*b) );
}


// =======================================================================================================
void calcStdDev( double *arr, double mean, double &std_dev )
{
    std_dev = 0;
    
    // stdDev:
    for ( int iSub = 0; iSub < nSubs; iSub++)
    {
        float diff = arr[iSub] - mean;
        std_dev += diff*diff;
    }
    
    // stdDev of the mean:
    std_dev /= (nSubs-1); // -1 is important. see also Marek's lecture - https://indico.cern.ch/event/722562/
    std_dev = sqrt(std_dev);
    std_dev = std_dev / sqrt(nSubs);
}


// =======================================================================================================
void      InitializeObjects();
void      /*ReadFitParamsFromTree*/ ReadFitParamsFromLineShapes(); //TString paramTreeName);
Double_t  EvalFitValue(Int_t particle, Double_t x);
void      RetrieveMoments(TIdentity2D *tidenObj, TVectorF *vecMom1st, TVectorF *vecMom2nd,  TVectorF *vecMom2ndMixed, TVectorF *vecInt);
void      PrintInitInfo();
//


// ======= Params =============================================================================
TString str_file_dEdx_in_bins = "file_dEdx_in_bins.root";
TString str_fitted_dEdx_in_pt_bins = "fitted_dEdx_in_pt_bins.root";


const Int_t fnEtaBins      = 1;//10;//1;
const Int_t fnCentBins     = 1;//9;
const Int_t fnParticleBins = 6;//4;//2;
const Int_t fnMomBins      = 130;

const Int_t fndEdxBins      = 2*300*5;//19;//120;
Int_t fMindEdx = -30000;
Int_t fMaxdEdx = 30000;
//
//
// Look up table related
const Int_t nBinsLineShape      = 400;
Int_t       fnTestEntries       = 0;
//
//
const Float_t fMomRangeDown = 0.2;//0.1;//0.3;
const Float_t fMomRangeUp   = 1.5;//2.0;//2.7;
const Int_t colors[]   = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed, kGreen};

//
// fixed tree branches --> [0]=event; [1]=dEdx; [2]=sign; [3]=cutBit; [4]=npart; [5]=Px; [6]=Py; [7]=Pz;
Double_t fTreeVariablesArray[10];//8];
const Int_t nBranches = 4;

TString branchNames[nBranches]={"centrality","pt", "eta", "mom"};
//
//
Int_t   fNthFitIteration = 6;
TString treeIdentity     = "treeIM";
//
//
// =======================================================================================================
//
// Inputs
Char_t  inputfileNameDataTree[255];     //file name of tree

TString fileNameDataTree = "";
Int_t fCentInputBin = 7;


//
Double_t nEvents = 0;
Double_t nnorm   = 1.;
Int_t fEtaBin, fCentBin, fMomBin;
Double_t fEtaForEval; // !!! IA spec for FB
Bool_t fFlagAfterLoop = false; // !!! IA spec for FB
Int_t fUsedBins[fnEtaBins][fnCentBins][fnMomBins];
UInt_t fCutBit;
TVectorF *fIntegrals;
TVectorF *fMoments1st,  *fMoments2nd, *fMoments2ndMixed;
//
// to be initialized
TH1D *fcRows   = NULL;
TH1D *fhPtot   = NULL;
TH1D *fhEta   = NULL;
static TH1D *****hParticles;
static TF1 *****fParticles;
TH1D **hDedxDebug;
TTree *momTree = NULL;
//
// members
TFile *fLineShapesLookUpTable = NULL;
TFile *fLineShapesLookUpTableFITTED = NULL;
TClonesArray *cloneArrHist=NULL;
TClonesArray *cloneArrFunc=NULL;
TTree *treeLookUp = NULL;
TFile *outFile = NULL;
TStopwatch timer;
//
Double_t fAmpArr[fnEtaBins][fnCentBins][fnMomBins][fnParticleBins]; // fAmpArr[fEtaBin][fCentBin][fMomBin][particleType]
Double_t fMeanArr[fnEtaBins][fnCentBins][fnMomBins][fnParticleBins];
Double_t fSigmaArr[fnEtaBins][fnCentBins][fnMomBins][fnParticleBins];
Double_t fSkewArr[fnEtaBins][fnCentBins][fnMomBins][fnParticleBins];
Double_t fKurtosisArr[fnEtaBins][fnCentBins][fnMomBins][fnParticleBins];


//
// -----------------------------------------------------------------------------------------
void analyse_identity_dEdx( int _cBin = 0 )
{
    sprintf(inputfileNameDataTree, "/Users/macbookpro/work/NICA_analysis/2020_07_SMASH_in_mpdroot/analysis_SMASH_output_2020_09_03/an_results_40_file_with_IM_tree_dEdx.root" ); // file with the identity tree
    fCentInputBin       = _cBin; //atof(argv[3]);
    // cout << fCentInputBin << endl;
    cout<<" main.Info: read file names from input "<<endl;
    TString outPutFileNAme = Form("TIMoments2D__cent_%d.root",fCentInputBin);
    outFile = new TFile(outPutFileNAme,"recreate");
    cout << " main.Info: write output into:    " << outPutFileNAme << endl;
    

    //
    InitializeObjects();
    fileNameDataTree   = inputfileNameDataTree;

    // Initialize objects and get the bin information
    TROOT IdentityMethod("IdentityMethod","compiled identity method");
    PrintInitInfo();
    ReadFitParamsFromLineShapes();//fileNameLineShapes);
    

    // Create the TIdentity2D object and start analysis
    TIdentity2D *iden4_all = new TIdentity2D(fnParticleBins);
    iden4_all -> SetBranchNames(nBranches,branchNames);
    iden4_all -> SetFileName(fileNameDataTree);
    iden4_all -> SetFunctionPointers(EvalFitValue);
    //    iden4_all -> SetLimits(fMindEdx,fMaxdEdx,250.,1500.,10); // --> (dEdxMin,dEdxMax,binwidth), if slice histograms are scaled wrt binwidth, then binwidth=1
    iden4_all -> SetLimits(fMindEdx,fMaxdEdx,fndEdxBins/*250.*/,1500.,10); // --> (dEdxMin,dEdxMax,binwidth), if slice histograms are scaled wrt binwidth, then binwidth=1
    //    iden4_all -> SetLimits(fMindEdx,fMaxdEdx,100.,150.,10); // --> (dEdxMin,dEdxMax,binwidth), if slice histograms are scaled wrt binwidth, then binwidth=1
    iden4_all -> SetUseSign(0);  // pass input sign value to TIdentity module
    iden4_all -> SetSeparateSign(kFALSE);
    Long_t nEntries;
    iden4_all -> GetTree(nEntries,treeIdentity);
    iden4_all -> Reset();
    
    
    TIdentity2D *arrIden4[nSubs];
    for( int subId = 0; subId < nSubs; subId++ )
    {
        arrIden4[subId] = new TIdentity2D(fnParticleBins);
        TIdentity2D *_iden = arrIden4[subId];
        
        _iden -> SetBranchNames(nBranches,branchNames);
        _iden -> SetFileName(fileNameDataTree);
        _iden -> SetFunctionPointers(EvalFitValue);
        //    _iden -> SetLimits(fMindEdx,fMaxdEdx,250.,1500.,10); // --> (dEdxMin,dEdxMax,binwidth), if slice histograms are scaled wrt binwidth, then binwidth=1
        _iden -> SetLimits(fMindEdx,fMaxdEdx,fndEdxBins/*250.*/,1500.,10); // --> (dEdxMin,dEdxMax,binwidth), if slice histograms are scaled wrt binwidth, then binwidth=1
        //    _iden -> SetLimits(fMindEdx,fMaxdEdx,100.,150.,10); // --> (dEdxMin,dEdxMax,binwidth), if slice histograms are scaled wrt binwidth, then binwidth=1
        _iden -> SetUseSign(0);  // pass input sign value to TIdentity module
        _iden -> SetSeparateSign(kFALSE);
        Long_t nEntries;
        _iden -> GetTree(nEntries,treeIdentity);
        _iden -> Reset();
    }
    //
    // track by track loop --> read all track info and add tracks to the iden4 object
    //
    
    if (fnTestEntries>0) nEntries = fnTestEntries;
    Int_t countEntry=0;
    int subsampleId = gRandom->Integer(nSubs);
    for( Int_t i = 0; i < nEntries; i++ )
    {
        
        // get entry for all
        bool haveEntry = iden4_all ->  GetEntry(i);
        
        if ( iden4_all->IsNextEvent() ) // NEXT EVENT! Change subsample number!
            subsampleId = gRandom->Integer(nSubs);
        
        if( /*!iden4_all ->  GetEntry(i)*/ !haveEntry ) continue;
        iden4_all      ->  GetBins(nBranches, fTreeVariablesArray);    // reads identity tree and retrives mybin[] info
        
        // get entry for subsample:
        TIdentity2D *_subs_iden4 = arrIden4[subsampleId];
        if( !_subs_iden4 ->  GetEntry(i) ) continue;
        _subs_iden4      ->  GetBins(nBranches, fTreeVariablesArray);    // reads identity tree and retrives mybin[] info
        
        //
        // Choose which kind of tree input is used
        fCutBit  = 0;
        fEtaBin  = 0;
        fCentBin = 0;
        Double_t pt = fTreeVariablesArray[5]; //sqrt(fTreeVariablesArray[5]*fTreeVariablesArray[5]+fTreeVariablesArray[6]*fTreeVariablesArray[6]);
        Double_t eta = fTreeVariablesArray[6];
        Double_t mom = fTreeVariablesArray[7];
        
        //        if(i%20000000 == 0) {
        if(i%100000 == 0) {
            cout << " main.Info: track " << i << " of " << nEntries;
            cout << " pt " << pt;
            cout << " -- bin0 =  " << fTreeVariablesArray[0];
            cout << " -- bin1 =  " << fTreeVariablesArray[1];
            cout << " -- bin2 =  " << fTreeVariablesArray[2];
            cout << " -- bin3 =  " << fTreeVariablesArray[3];
            cout << " -- bin4 =  " << fTreeVariablesArray[4] << endl;
        }
                
        // Forward-Backward correlations analysis in eta-windows (−0.5, −0.1) − (0.1, 0.5):
        if( fabs(eta) > 0.5 || fabs(eta) < 0.1 ) // for Forward-Backward correlations analysis
            continue;
        
        // mom cuts
        if( mom < 0.3 || mom > 1.5 )//1.2)//2)
            continue;
        
        fEtaForEval = eta;
        
        //        cout << " -- dca =  " << dca << endl;
        
        fCentBin = 0; //fTreeVariablesArray[4];
        fMomBin  = fhPtot -> FindBin(/*pt*/mom) -1;
        
        // eta bin:
        fEtaBin  = 0; //fhEta -> FindBin(eta) -1;
        //        cout << fEtaBin << endl;
        
        
        fcRows->Fill(pt);

        // cout << "momentum " << momentum << " " << fMomBin << " " << fTreeVariablesArray[1] << endl;
        // cout << "fCentBin " << fCentBin << " " << fCentInputBin << endl;

        if( fCentBin != fCentInputBin )
            continue;
        if(fMomBin < 0 || fMomBin > fnMomBins -1)
        {
            cout << "AHTUNG! fMomBin = " << fMomBin << endl;
            continue;
        }
        //        hDedxDebug[fMomBin]->Fill( dca/*dca*/ );
        //
        // cout << "fCentBin " << fCentBin << " fMomBin " << fMomBin << endl;
        fUsedBins[fEtaBin][fCentBin][fMomBin] = 1;

        iden4_all -> AddEntry();
        _subs_iden4 -> AddEntry();

        countEntry++;
    } // end of track loop
    cout << "TEST after the loop" << endl;
    
    fFlagAfterLoop = true;

    //
    // Dump some debug histograms
    outFile->cd();
    for (Int_t i=0; i<fnMomBins/*120*/; i++)
    {
        hDedxDebug[i]->Write();
        for (Int_t j=0; j<fnParticleBins; j++) {
            // cout << "fCentBin " << fEtaBin << " " << fCentBin << " " << i << " " << j << endl;
            hParticles[fEtaBin][fCentInputBin][i][j]->Write();
        }
    }
    //
    //
    cout << "main.Info: Total number of tracks processed = " << countEntry << endl;
    iden4_all -> Finalize();
    for( int subId = 0; subId < nSubs; subId++ )
        arrIden4[subId] -> Finalize();
    
    //
    // Calculate 2. order moments only for full range
    cout << " ==================================" << endl;
    cout << " main.Info: calculating integrals " <<endl;
    timer.Reset(); timer.Start();
    cout << " ==================================" << endl;
    for(Int_t i = 0; i < fnEtaBins; i++){
        for(Int_t j = 0; j < fnCentBins; j++){
            for(Int_t k = 0; k < fnMomBins; k++){
                
                if(fUsedBins[i][j][k] != 1) continue;
                fEtaBin   =  i;
                fCentBin  =  j;
                fMomBin   =  k;
                iden4_all  -> AddIntegrals(0); // real sign information passed for the check with real data tree
                for( int subId = 0; subId < nSubs; subId++ )
                    arrIden4[subId] -> AddIntegrals(0);
            }
        }
    }
    iden4_all -> CalcMoments();
    for( int subId = 0; subId < nSubs; subId++ )
        arrIden4[subId] -> CalcMoments();
    //
    // Retrive Moments
    RetrieveMoments(iden4_all,fMoments1st,fMoments2nd,fMoments2ndMixed,fIntegrals);
    
    // subsamples:
    const int nSpecies = fnParticleBins;
    double _subsN[nSpecies][nSubs] = {};
    double _subsN2[nSpecies][nSubs] = {};
    double _subsCross[nSpecies*(nSpecies-1)/2][nSubs] = {};
    double _subs_bcorr[nSpecies*(nSpecies-1)/2][nSubs] = {};
    
    for( int subId = 0; subId < nSubs; subId++ )
    {
        cout << "##### retrieving moments for subId = " << subId << endl;
        //        RetrieveMoments(arrIden4[subId],fMoments1st,fMoments2nd,fMoments2ndMixed,fIntegrals);
        
        //        TVectorF *vecMom1st,  *vecMom2nd, *vecMom2ndMixed;
        
        TIdentity2D *tidenObj = arrIden4[subId];
        // 1st Moments
        for (Int_t i=0; i<fnParticleBins; i++){
            double mom = tidenObj -> GetMean(i);
            //            cout << " First  --> " << i << " = " <<  mom << endl;
            //            histFirstMom->SetBinContent(i+1,  mom);
            _subsN[i][subId] = mom;
            //        histFirstMom->SetBinError(i+1, std_dev_N[i]);
        }
        // 2nd Moments
        for (Int_t i=0; i<fnParticleBins; i++){
            double mom = tidenObj -> GetSecondMoment(i);
            //            cout << " Second --> " << i << " = " <<  mom << endl;
            //            histSecondMom->SetBinContent(i+1, mom;
            _subsN2[i][subId] = mom;
            //        histSecondMom->SetBinError(i+1, std_dev_N2[i]);
        }
        // Mixed Moments
        int counterCross = 0;
        for (Int_t i=0; i<fnParticleBins; i++){
            for (Int_t j=0; j<fnParticleBins; j++){
                if (i>=j) continue;
                double mom = tidenObj -> GetMixedMoment(i,j);
                //                cout << " Mixed --> " << i << "-" << j << " = " <<  mom << endl;
                //                histMixedMom->SetBinContent(counterCross+1, mom );
                _subsCross[counterCross][subId] = mom;
                
                _subs_bcorr[counterCross][subId] =  calc_bcorr( mom, _subsN[i][subId], _subsN[j][subId], _subsN2[i][subId], _subsN2[j][subId] );
                
                //            histMixedMom->SetBinError(objCounter+1, std_dev_Cross[objCounter] );
                counterCross++;
            }
        }
    } // end of loop over subsamples
    
    double av_bcorr[nSpecies*(nSpecies-1)/2] = {};
    
    // calc std dev:
    double std_dev_N[nSpecies] = {};
    double std_dev_N2[nSpecies] = {};
    double std_dev_Cross[nSpecies*(nSpecies-1)/2] = {};
    double std_dev_bcorr[nSpecies*(nSpecies-1)/2] = {};
    
    int counterCross = 0;
    for( int i = 0; i < nSpecies; i++ )
    {
        calcStdDev( _subsN[i], /*avN[i]*/(*fMoments1st)[i], std_dev_N[i] );
        calcStdDev( _subsN2[i], /*avN2[i]*/(*fMoments2nd)[i], std_dev_N2[i] );
        for( int j = i+1; j < nSpecies; j++ )
        {
            calcStdDev( _subsCross[counterCross], /*avCross[counterCross]*/(*fMoments2ndMixed)[counterCross], std_dev_Cross[counterCross] );
            
            av_bcorr[counterCross] = calc_bcorr( (*fMoments2ndMixed)[counterCross],  (*fMoments1st)[i], (*fMoments1st)[j], (*fMoments2nd)[i], (*fMoments2nd)[j] );
            
            calcStdDev( _subs_bcorr[counterCross], av_bcorr[counterCross], std_dev_bcorr[counterCross] );
            
            counterCross++;
        }
    }
    
    
    // ### Prepare histos for reconstructed moments:
    
    // ### histos with moments
    //    const int nSpecies = fnParticleBins;
    TH1D *histFirstMom = new TH1D( "histFirstMom","; ;Entries",nSpecies,-0.5,nSpecies-0.5);
    TString arrFirstMomNames[nSpecies] = {"#LT #pi^{+} #GT", "#LT K^{+} #GT", "#LT p #GT", "#LT #pi^{#minus} #GT", "#LT K^{#minus} #GT", "#LT #bar{p} #GT" };
    for(Int_t i = 0; i < nSpecies; i++)
        histFirstMom->GetXaxis()->SetBinLabel( i+1, arrFirstMomNames[i].Data() );
    //    tuneHist1D(histFirstMom);
    
    TH1D *histSecondMom = new TH1D( "histSecondMom","; ;Entries", nSpecies, -0.5, nSpecies-0.5 );
    TString arrSecMomNames[nSpecies] = {"#LT #pi^{+ 2} #GT", "#LT K^{+ 2} #GT", "#LT p^{2} #GT", "#LT #pi^{#minus 2} #GT", "#LT K^{#minus 2} #GT", "#LT #bar{p}^{2} #GT" };
    for(Int_t i = 0; i < nSpecies; i++)
        histSecondMom->GetXaxis()->SetBinLabel( i+1, arrSecMomNames[i].Data() );
    //    tuneHist1D(histSecondMom);
    
    const int nComb = nSpecies*(nSpecies-1)/2;
    TH1D *histMixedMom = new TH1D( "histMixedMom","; ;Entries", nComb, -0.5, nComb-0.5 );
    TString arrMixedMomNames[nComb] = { "#LT #pi^{+}K^{+} #GT", "#LT #pi^{+}p #GT", "#LT #pi^{+}#pi^{#minus} #GT", "#LT #pi^{+}K^{#minus} #GT", "#LT #pi^{+}#bar{p} #GT",
                                        "#LT K^{+}p #GT",   "#LT K^{+}#pi^{#minus} #GT", "#LT K^{+}K^{#minus} #GT", "#LT K^{+}#bar{p} #GT",
                                        "#LT p#pi^{#minus} #GT", "#LT pK^{#minus} #GT", "#LT p#bar{p} #GT",
                                        "#LT #pi^{#minus}K^{#minus} #GT", "#LT #pi^{#minus}#bar{p} #GT",
                                        "#LT K^{#minus}#bar{p} #GT",
                                      };
    for(Int_t i = 0; i < nComb; i++)
        histMixedMom->GetXaxis()->SetBinLabel( i+1, arrMixedMomNames[i].Data() );
    //    tuneHist1D(histMixedMom);
    
    
    // bcorr:
    TH1D *hist_bcorr = new TH1D( "hist_bcorr","; ;Entries", nComb, -0.5, nComb-0.5 );
    TString arrBcorrBinsNames[nComb] = { "#pi^{+}K^{+}", "#pi^{+}p", "#pi^{+}#pi^{#minus}", "#pi^{+}K^{#minus}", "#pi^{+}#bar{p}",
                                         "K^{+}p",   "K^{+}#pi^{#minus}", "K^{+}K^{#minus}", "K^{+}#bar{p}",
                                         "p#pi^{#minus}", "pK^{#minus}", "p#bar{p}",
                                         "#pi^{#minus}K^{#minus}", "#pi^{#minus}#bar{p}",
                                         "K^{#minus}#bar{p}",
                                       };
    for(Int_t i = 0; i < nComb; i++)
        hist_bcorr->GetXaxis()->SetBinLabel( i+1, arrBcorrBinsNames[i].Data() );
    //    tuneHist1D(hist_bcorr);
    
    
    
    histFirstMom->SetMarkerStyle(25);
    histFirstMom->SetMarkerColor(kRed);
    histFirstMom->SetLineColor(kRed);
    
    histSecondMom->SetMarkerStyle(25);
    histSecondMom->SetMarkerColor(kRed);
    histSecondMom->SetLineColor(kRed);
    
    histMixedMom->SetMarkerStyle(25);
    histMixedMom->SetMarkerColor(kRed);
    histMixedMom->SetLineColor(kRed);
    
    
    hist_bcorr->SetMarkerStyle(25);
    hist_bcorr->SetMarkerColor(kRed);
    hist_bcorr->SetLineColor(kRed);
    
    
    
    
    // #################  PRINTING WITH ERRORS:
    cout << "##### PRINTING WITH ERRORS: " << endl;
    // 1st Moments
    for (Int_t i=0; i<fnParticleBins; i++)
    {
        cout << " First  --> " << i << " = " <<  (*fMoments1st)[i] << " +/- " <<  std_dev_N[i] << endl;
        histFirstMom->SetBinContent(i+1,  (*fMoments1st)[i]);
        histFirstMom->SetBinError(i+1, std_dev_N[i]);
    }
    // 2nd Moments
    for (Int_t i=0; i<fnParticleBins; i++){
        //        (*vecMom2nd)[i] = tidenObj -> GetSecondMoment(i);
        cout << " Second --> " << i << " = " <<  (*fMoments2nd)[i] << " +/- " <<  std_dev_N2[i] << endl;
        histSecondMom->SetBinContent(i+1, (*fMoments2nd)[i]);
        histSecondMom->SetBinError(i+1, std_dev_N2[i]);
    }
    // Mixed Moments
    //    Int_t objCounter = 0;
    counterCross = 0;
    for (Int_t i=0; i<fnParticleBins; i++){
        for (Int_t j=0; j<fnParticleBins; j++){
            if (i>=j) continue;
            cout << " Mixed --> " << i << "-" << j << " = " <<  (*fMoments2ndMixed)[counterCross] << " +/- " <<  std_dev_Cross[counterCross] << endl;
            histMixedMom->SetBinContent(counterCross+1, (*fMoments2ndMixed)[counterCross] );
            histMixedMom->SetBinError(counterCross+1, std_dev_Cross[counterCross] );
            counterCross++;
        }
    }
    
    
    counterCross = 0;
    for (Int_t i=0; i<fnParticleBins; i++){
        for (Int_t j=0; j<fnParticleBins; j++){
            if (i>=j) continue;
            cout << " b_corr --> " << i << "-" << j << " = " <<  av_bcorr[counterCross] << " +/- " <<  std_dev_bcorr[counterCross] << endl;
            hist_bcorr->SetBinContent(counterCross+1, av_bcorr[counterCross] );
            hist_bcorr->SetBinError(counterCross+1, std_dev_bcorr[counterCross] );
            counterCross++;
        }
    }
    
    
    
    
    TCanvas *canv_FirstMom = new TCanvas("canv_FirstMom","FirstMom",120,80,800,600 );
    //    tuneCanvas( canv_FirstMom );
    histFirstMom->DrawCopy();
    
    
    TCanvas *canv_SecMom = new TCanvas("canv_SecMom","SecMom",150,110,800,600 );
    //    tuneCanvas( canv_SecMom );
    histSecondMom->DrawCopy();
    
    
    TCanvas *canv_MixedMom = new TCanvas("canv_MixedMom","MixedMom",170,120,800,600 );
    //    tuneCanvas( canv_MixedMom );
    histMixedMom->DrawCopy();
    
    TCanvas *canv_bcorr = new TCanvas("canv_bcorr","bcorr",190,140,800,600 );
    //    tuneCanvas( canv_bcorr );
    hist_bcorr->DrawCopy();
    
    // ### write results (moments)
    TFile *file_extracted_moments_IM = new TFile( "file_extracted_moments_IM.root", "recreate" );
    histFirstMom->Write();
    histSecondMom->Write();
    histMixedMom->Write();
    hist_bcorr->Write();
    file_extracted_moments_IM->Close();
    
    
    
    //
    cout << "====================================" << endl;
    cout << " main.Info: calculation is finished " << endl;
    timer.Stop(); timer.Print();
    cout << "====================================" << endl;
    //
    // Fill output tree
    cout << " main.Info: fill tree " << endl;
    momTree -> Fill();
    //
    // Close file and clear memory
    cout << " main.Info: dump output " << endl;
    outFile -> cd();
    fcRows  -> Write("fcRows");
    momTree -> Write();
    outFile -> Close();
    delete outFile; //yeni eklave etdim.
    delete iden4_all;
    return 1;
}


// -----------------------------------------------------------------------------------------
void ReadFitParamsFromLineShapes()
{
    // Read the lineShapes from the lookup table and fill them into pointer arrays
    
    cout << " --- In ReadFitParamsFromLineShapes --- " << endl;
    fLineShapesLookUpTable = new TFile( str_file_dEdx_in_bins );
    fLineShapesLookUpTableFITTED = new TFile( str_fitted_dEdx_in_pt_bins );
    
    for (Int_t eBin = 0; eBin<fnEtaBins; eBin++)
        for (Int_t cBin = 0; cBin<fnCentBins; cBin++)
        {
            for (Int_t pid = 0; pid<fnParticleBins; pid++)
            {
                TH2D *h2D_Particles = (TH2D*) fLineShapesLookUpTable->Get( Form( "fDeDx_%d", pid ) );
                TH1D *proj = h2D_Particles->ProjectionX();
                for (Int_t ptBin = 0; ptBin<fnMomBins; ptBin++)
                {
                    // cout << eBin << " " << cBin << " " << ptBin << " " << pid << endl;
                    double lowEdge = proj->GetBinLowEdge(ptBin+1);
                    double upEdge = lowEdge + proj->GetBinWidth(ptBin+1);
                    
                    
                    hParticles[eBin][cBin][ptBin][pid] = (TH1D*) fLineShapesLookUpTableFITTED->Get( Form( "fDeDx_eta%d_centr%d_pt%d_range_%.2f_%.2f_pid%d", eBin, cBin, ptBin, lowEdge, upEdge, pid ) );
                    //cout << hParticles[eBin][cBin][ptBin][pid] << endl;
                    
                    TH1D *h = hParticles[eBin][cBin][ptBin][pid];
                    TF1 *func = h->GetFunction("fGenGaus");
                    
                    if(func)
                    {
                        for ( int i = 0; i < h->GetNbinsX(); i++ )
                        {
                            double value = func->Eval( h->GetBinCenter(i+1) );
                            h->SetBinContent( i+1, value );
                        }
                        
                    }
                }
            }
        }
}


// -----------------------------------------------------------------------------------------
Double_t EvalFitValue(Int_t particle, Double_t x)
{
    if( !fFlagAfterLoop )
    {
        int winFB = (fEtaForEval>0) ? +1 : -1; //( fEtaBin < fnEtaBins/2 ) ? -1 : +1; // F or B
        //    cout << "EvalFitValue" << endl;
        if ( particle < 3  /*F*/ && winFB == -1 )   return 0;   // because particle is in B window
        if ( particle >= 3 /*B*/ && winFB == 1  )   return 0;   // because particle is in F window
    }
        
    Int_t bin = hParticles[fEtaBin][fCentBin][fMomBin][particle]->FindBin(x);
    //    cout << "bin = " << bin << endl;
    // cout << fEtaBin << " " << fCentBin <<  " " << fMomBin <<  " " << particle <<  " " << x << " " << fParticles[fEtaBin][fCentBin][fMomBin][particle]->Eval(x) << endl;
    //    cout << fEtaBin << " " << fCentBin <<  " " << fMomBin <<  " pid=" << particle <<  " dEdx=" << x << " bin=" << bin << " content=" << hParticles[fEtaBin][fCentBin][fMomBin][particle]->GetBinContent(bin) << endl;
    return hParticles[fEtaBin][fCentBin][fMomBin][particle]->GetBinContent(bin);
}




// -----------------------------------------------------------------------------------------
void InitializeObjects()
{
    //
    fcRows  = new TH1D("fcRows","fcRows",240 ,0, 2.4 );
    fhPtot  = new TH1D("fhPtot","Momentum Bins"  ,fnMomBins ,fMomRangeDown, fMomRangeUp );
    fhEta  = new TH1D("fhEta","Eta Bins"  , 10, -1.0, 1.0 );
    //
    hDedxDebug = new TH1D *[fnMomBins];
    for (Int_t i = 0; i < fnMomBins; i++) {
        hDedxDebug[i] = new TH1D(Form("hDedx_%d", i),Form("hDedx_%d", i),nBinsLineShape,fMindEdx,fMaxdEdx);
    }
    //
    fParticles = new TF1 ****[fnEtaBins];
    for (Int_t i = 0; i<fnEtaBins; i++)
    {
        fParticles[i] = new TF1***[fnCentBins];
        for (Int_t j = 0; j<fnCentBins; j++){
            fParticles[i][j] = new TF1**[fnMomBins];
            for (Int_t k = 0; k<fnMomBins; k++){
                fParticles[i][j][k] = new TF1*[fnParticleBins];
                for (Int_t ii = 0; ii<fnParticleBins; ii++){
                    fParticles[i][j][k][ii] = NULL;
                }
            }
        }
    }
    hParticles = new TH1D ****[fnEtaBins];
    for (Int_t i = 0; i<fnEtaBins; i++)
    {
        hParticles[i] = new TH1D***[fnCentBins];
        for (Int_t j = 0; j<fnCentBins; j++){
            hParticles[i][j] = new TH1D**[fnMomBins];
            for (Int_t k = 0; k<fnMomBins; k++){
                hParticles[i][j][k] = new TH1D*[fnParticleBins];
                for (Int_t ii = 0; ii<fnParticleBins; ii++){
                    hParticles[i][j][k][ii] = NULL;
                }
            }
        }
    }
    //
    //
    for( Int_t i = 0; i < fnEtaBins; i++ )
    {
        for( Int_t j = 0; j < fnCentBins; j++ ){
            for( Int_t k = 0; k < fnMomBins; k++ ){
                fUsedBins[i][j][k] = -1;
            }
        }
    }
    //
    // initialize counters
    fIntegrals    = new TVectorF(fnParticleBins);
    fMoments1st   = new TVectorF(fnParticleBins);
    fMoments2nd   = new TVectorF(fnParticleBins);
    fMoments2ndMixed   = new TVectorF(fnParticleBins*fnParticleBins);
    
    for(Int_t i=0;i<fnParticleBins*fnParticleBins; i++)
    {
        if (i<fnParticleBins)
        {
            (*fIntegrals)[i]=0.;
            (*fMoments1st)[i]=0.;
            (*fMoments2nd)[i]=0.;
        }
        (*fMoments2ndMixed)[i]=0.;
    }
    //
    // initialize output tree
    //
    momTree = new TTree("momTree","momTree");
    momTree -> Branch("nEvents",&nEvents);
    momTree -> Branch("nnorm",&nnorm);
    momTree -> Branch("fIntegrals",&fIntegrals);
    momTree -> Branch("fMoments1st",&fMoments1st);
    momTree -> Branch("fMoments2nd",&fMoments2nd);
    momTree -> Branch("fMoments2ndMixed",&fMoments2ndMixed);
    
}


// -----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------
void RetrieveMoments(TIdentity2D *tidenObj, TVectorF *vecMom1st, TVectorF *vecMom2nd,  TVectorF *vecMom2ndMixed, TVectorF *vecInt)
{
    
    // 1st Moments
    for (Int_t i=0; i<fnParticleBins; i++){
        (*vecMom1st)[i] = tidenObj -> GetMean(i);
        cout << " First  --> " << i << " = " <<  (*vecMom1st)[i] << endl;
        //        histFirstMom->SetBinContent(i+1,  (*vecMom1st)[i]);
        //        histFirstMom->SetBinError(i+1, std_dev_N[i]);
    }
    // 2nd Moments
    for (Int_t i=0; i<fnParticleBins; i++){
        (*vecMom2nd)[i] = tidenObj -> GetSecondMoment(i);
        cout << " Second --> " << i << " = " <<  (*vecMom2nd)[i] << endl;
        //        histSecondMom->SetBinContent(i+1, (*vecMom2nd)[i]);
        //        histSecondMom->SetBinError(i+1, std_dev_N2[i]);
    }
    // Mixed Moments
    //    Int_t objCounter = 0;
    int counterCross = 0;
    for (Int_t i=0; i<fnParticleBins; i++){
        for (Int_t j=0; j<fnParticleBins; j++){
            if (i>=j) continue;
            (*vecMom2ndMixed)[counterCross] = tidenObj -> GetMixedMoment(i,j);
            cout << " Mixed --> " << i << "-" << j << " = " <<  (*vecMom2ndMixed)[counterCross] << endl;
            //            histMixedMom->SetBinContent(counterCross+1, (*vecMom2ndMixed)[counterCross] );
            //            histMixedMom->SetBinError(objCounter+1, std_dev_Cross[objCounter] );
            counterCross++;
        }
    }
    
    if( fnParticleBins>2 )
    {
        cout << "#########" << endl;
        
        for (Int_t i=0; i<fnParticleBins; i++)
        {
            double Var = (*vecMom2nd)[i] - (*vecMom1st)[i]*(*vecMom1st)[i];
            cout << "PID " << i << ": Var = " << Var << endl;
        }
        
        for (Int_t i=0; i<fnParticleBins; i++)
        {
            for (Int_t j=i+1; j<fnParticleBins; j++)
            {
                double VarA = (*vecMom2nd)[i] - (*vecMom1st)[i]*(*vecMom1st)[i];
                double VarB = (*vecMom2nd)[j] - (*vecMom1st)[j]*(*vecMom1st)[j];
                double COV_FB = tidenObj -> GetMixedMoment(i,j) - (*vecMom1st)[i] * (*vecMom1st)[j];
                double bcorr = COV_FB / sqrt(VarA*VarB);
                cout << "COV_FB (" << i << ", " << j << ") = " << COV_FB << ", bcorr = " << bcorr << endl;
            }
            
        }
    }
    
    //
    //Integrals:
    // 1st Moments
    for (Int_t i=0; i<fnParticleBins; i++){
        (*vecInt)[i] = tidenObj -> GetMeanI(i);
        cout << " Integrals --> " << i << " = " <<  (*vecInt)[i] << endl;
    }
}


// -----------------------------------------------------------------------------------------
void PrintInitInfo()
{
    
    //
    cout << " ================================================================================= " << endl;
    cout << " InitializeObjects.Info: treeIdentity          = " << treeIdentity       << endl;
    cout << " ================================================================================= " << endl;
    cout << " InitializeObjects.Info: Inputs: " << endl;
    cout << " InitializeObjects.Info: data Tree             = " << fileNameDataTree       << endl;
    cout << " ================================================================================= " << endl;
    //
}
