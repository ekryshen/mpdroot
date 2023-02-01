
#include "TpcSectorGeoAlignmentVAK.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "RecoMiniMC.h"
#include "miniTpcAlignment.h"
#include "TpcAlignment.h"
#include <Riostream.h>
#include <cstdio>                               // для функции remove

//common variables
Double_t R0_A[3][24],a_A[24],b_A[24],g_A[24];
Long64_t NhitSec0[24];
Double_t SecChi2_0[24];
Int_t maxS,ns=0;
Double_t kdAr=1.,kdAa=2.;
Double_t kdAR=kdAr,kdAA=kdAa;
Double_t gdX[3];

void getAlignment(const char* afile);
void getRealAlignment(const char* afile);
void putAlignment(const char* afile);
void convert2miniMC(const char* InDataFile,const char* InAlignmentFile);
inline void cleanfiles(const char* OutFile,const char* AFile) {
  printf("rm OutFile is %d    rm AFile is %d\n",remove(OutFile),remove(AFile));
}
Double_t getChi2(const char* infile);
Double_t getNhits(const char* infile);
Double_t checkdiv(Double_t* dx,Double_t* dy,Double_t* dz,
                  Double_t* da,Double_t* db,Double_t* dg,Int_t* Sectors);
//  Int_t Sectors[24]={1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0};
//                   0         5        10        15        20
//                   |         |         |         |         |
  Int_t Sectors[24]={1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0};

// = = = = = = = = = = = = = 

int tpcAlignment(const char* InDataFile,const char* InAlignmentFile=" "){

// macro to find the TPC alignment using experimental track hits
// 1.Input file (inDataFile) with track hits and the information about 
//   magnetic and electric filds
// 2.Input file (InAlignmentFile) with the alignment used at the reconstruction.
//   If the name is " ", then file inDataFile has the miniMC structure
//   and the real detector alignment is known.
  Int_t nLoop=1000;                                // maximum itteration number
  Int_t evt1=0,  evt2=10000;                      // first & last tracks
  Double_t Epsilon=1e-5;                           // delta F limit 
  Int_t recoPL=0,  alignmentPL=0;                  // printing levels
  Int_t retN=3;                                    // attemps to leave local
                                                   // minmum  
  Int_t MinHits = 20;                              // minimum hits number
  Double_t EvalCut=4;                              // Minuit2 maximum EvalCut 
  
//  input data
  const char* InFile; 
  TString inData=InDataFile;
  TString inAlignment=InAlignmentFile;
  if(!inAlignment.EqualTo(" ")) {        // real DST data format
    // convert2miniMC(InDataFile,InAlignmentFile);
    printf("Do not know about real DST format for the Alignment task, stop");
    return -2;
  }
  else {
    InFile=InDataFile; 
  }

  const char* oDir= " ";                            //  output dir for Reco
     
  TString tmp="temp";
  TString inF=InFile;
  Ssiz_t lif=inF.Length();
  Ssiz_t ld=inF.Last('/');
  stringstream stempOut;
  stringstream stoutA;
  if(ld>=0) {
    TString prename=inF(0,ld+1);         // variant name = aFile-DirName
    TString ownname=inF(ld+1,lif-1);         // variant name = aFile-DirName
    stempOut<<prename.Data()<<tmp.Data()<<"_"<<ownname.Data();
    stoutA<<tmp.Data()<<"A_"<<ownname.Data();
//    stempOut<<prename.Data()<<"temp_"<<ownname.Data();
//    stoutA<<"tempA_"<<ownname.Data();
  } else {
    stempOut<<tmp.Data()<<"_"<<InFile;
    stoutA<<tmp.Data()<<"A_"<<InFile;
  }
  TString* nfile=new TString(stempOut.str());
  const char* OutFile=nfile->Data();
  TString* afile=new TString(stoutA.str());
  const char* AFile=afile->Data();               // output alignment file name
 
//  Int_t Sectors[24]={1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0};
//                   0         5        10        15        20
//                   |         |         |         |         |
  Int_t Sectors[24]={1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0};
  Long64_t NhitSec[24],Nhits;
  Double_t SecChi2[24];
  Double_t F0T=0,F0,F,sum,sum0,sum_0,eps0,eps;
  const char* Comment="";

  Double_t r0[3][24],a0[24],b0[24],g0[24];          // initial alignment values
  Double_t rr0[3][24],ra[24],rb[24],rg[24];          // real alignment values
  Double_t dRx[24],dRy[24],dRz[24],da[24],db[24],dg[24];
  Double_t dR0x[24],dR0y[24],dR0z[24],da0[24],db0[24],dg0[24];
  Double_t R0[3][24],dR0[3][24];
  Double_t sX[3]={0,0,0},sA[3]={0,0,0};

  BaseTpcSectorGeo *geosec = new TpcSectorGeoAlignmentVAK();   // Zero alignment
//  geosec->PrintSectorAlignment(0);
//  miniTpcAlignment* al=new miniTpcAlignment();
  TpcAlignment* al=new TpcAlignment();
  Int_t AbsMinHits=TMath::Abs(MinHits);
      RecoMiniMC* reco=new RecoMiniMC();

  for(Int_t s=0;s<24;s++) if(Sectors[s]) ns++;
  {                          // start of the real alignment search

    cout<<"==========================  files =========================="<<endl;
    cout<<"input file:           "<<InFile<<endl;
    cout<<"output reco file:     "<<OutFile<<endl;
    cout<<"output alignment file:  "<<AFile<<endl;
    cout<<"============================================================"<<endl;
    cout<<endl;

    cout<<"     maximum number of loops is "<<nLoop<<"  eps="<<Epsilon
      <<"  events="<<evt2-evt1<<endl;
    cout<<endl;
    if(inAlignment.EqualTo(" ")) {        // miniMC data
      // get real alignment we have to find and put it into <R0_A,a_A,b_A,g_A>
      getRealAlignment(InFile);
      for(Int_t s=0;s<24;s++) {
        if(Sectors[s]) {
          for(Int_t k=0;k<3;k++) rr0[k][s]=R0_A[k][s];
          ra[s]=a_A[s];rb[s]=b_A[s];rg[s]=g_A[s];
        }
      }
      putAlignment(AFile);
///      RecoMiniMC* reco=new RecoMiniMC();
      reco->SetPrintLevel(0);
      reco->Init(*geosec, evt1,evt2,Sectors,AbsMinHits,InFile,AFile,oDir,Comment,OutFile);
      reco->SetEvalCut(2*EvalCut);
      reco->reco();
///      delete reco;
      F0T=getChi2(OutFile);
    }

    // initial F0 for needed number of events 
    getAlignment(InFile);
    if(evt2!=0) {
      putAlignment(AFile);
///      RecoMiniMC* reco=new RecoMiniMC();
      reco->SetPrintLevel(recoPL);
      reco->Init(*geosec, evt1,evt2,Sectors,AbsMinHits,InFile,AFile,oDir,Comment,OutFile);
      reco->SetEvalCut(EvalCut);
      reco->reco();
///      delete reco;
      F0=getChi2(OutFile);
      Nhits=getNhits(OutFile);
    } else {
      F0=getChi2(InFile);
      Nhits=getNhits(InFile);
    }
    if(F0T==0) F0T=F0;
    if(inAlignment.EqualTo(" "))         // miniMC data
      printf("\n<<<<<<<<<<<<<<<<< F0T=%17.15f >>>>>>>>>>>>>>>>>>\n\n",F0T);
    // get initial alignment
    printf("========== initial data ==== events(%d-%d) ===== F0=%17.15f ====\n",
           evt1,evt2,F0);
    for(Int_t s=0;s<24;s++) {
      if(Sectors[s]) {
        for(Int_t k=0;k<3;k++) r0[k][s]=R0_A[k][s];
        a0[s]=a_A[s];b0[s]=b_A[s];g0[s]=g_A[s];
        if(inAlignment.EqualTo(" ")) {        // miniMC data
          printf("to find: sector=%2d  R0(%11.4e,%11.4e,%11.4e)  ",
                 s,rr0[0][s],rr0[1][s],rr0[2][s]);
          printf("abg(%11.4e,%11.4e,%11.4e)\n",ra[s],rb[s],rg[s]);
        }
        printf("initial: sector=%2d  R0(%11.4e,%11.4e,%11.4e)  ",
               s,R0_A[0][s],R0_A[1][s],R0_A[2][s]);
        printf("abg(%11.4e,%11.4e,%11.4e)\n",a_A[s],b_A[s],g_A[s]);
       for(Int_t k=0;k<3;k++) sX[k]+=rr0[k][s];
       sA[0]+=ra[s];sA[2]+=rb[s];sA[2]+=rg[s];
      }
    }
    printf("  initial mean sifts(%11.4e,%11.4e,%11.4e  %11.4e,%11.4e,%11.4e)\n",
           sX[0]/ns,sX[1]/ns,sX[2]/ns,sA[0]/ns,sA[1]/ns,sA[2]/ns);
    printf("===++++= kdAr=%f kdAa=%f ========== Nevt=%d ===========\n\n",
           kdAr,kdAa,evt2-evt1);
    for(Int_t k=0;k<3;k++) sX[k]=0;
    sA[0]=0;sA[2]=0;sA[2]=0;

    // find the step to the minimum
    al->Init(*geosec,evt1,evt2,AbsMinHits,Sectors,InFile," ");
    al->SetPrintLevel(alignmentPL);
    al->alignment();
    al->results(dRx,dRy,dRz,da,db,dg);
    checkdiv(dRx,dRy,dRz,da,db,dg,Sectors);
///    delete al;
  }

  // -> the best alignment loop
  // ******************************
  Int_t retC=retN;
  for(Int_t loop=1;loop<=nLoop;loop++) {
    cout<<"============== loop="<<loop<<"===========\n";
    for(Int_t s=0;s<24;s++) 
      if(Sectors[s]) {
        R0_A[0][s]+=kdAr*dRx[s];R0_A[1][s]+=kdAr*dRy[s];R0_A[2][s]+=kdAr*dRz[s];
        a_A[s]+=kdAa*da[s];b_A[s]+=kdAa*db[s];g_A[s]+=kdAa*dg[s];
        printf("sector=%2d  dR(%11.4e,%11.4e,%11.4e)  dA(%11.4e,%11.4e,%11.4e)\n",
               s,dRx[s],dRy[s],dRz[s],da[s],db[s],dg[s]);
        for(Int_t k=0;k<3;k++) sX[k]+=kdAr*gdX[k];
        sA[0]+=kdAa*da[s];sA[1]+=kdAa*db[s];sA[2]+=kdAa*dg[s];
      }
    putAlignment(AFile);

    // calculate Chi2 for the corrected alignment
///    RecoMiniMC* reco=new RecoMiniMC();
    reco->SetPrintLevel(recoPL);
    reco->Init(*geosec, evt1,evt2,Sectors,MinHits,InFile,AFile,oDir,Comment,OutFile);
    reco->SetEvalCut(2*EvalCut);
    reco->reco();
///    delete reco;
    F=getChi2(OutFile);
    Nhits=getNhits(OutFile);
    eps=TMath::Abs((F0-F)/F0);
    // print intermediate results
    printf("          chi2/ndf sum %7.5f  nhits=%lld\n\n",F,Nhits);
    if(inAlignment.EqualTo(" "))         // miniMC data
      printf("******** F(F0)=%17.15f(%17.15f)  eps(epsT)=%9.3e(%9.3e) ********\n",
            F,F0,eps,(F-F0T)/F0T);
    else
      printf("******** F(F0)=%17.15f(%17.15f)  eps=%9.3e ********\n",
            F,F0,eps);
    for(Int_t s=0;s<12;s++) {
      if(Sectors[s]) {
        printf("sector=%2d",s);
        printf("  R0(%11.4e,%11.4e,%11.4e)  abg(%11.4e,%11.4e,%11.4e)\n",
               R0_A[0][s],R0_A[1][s],R0_A[2][s],a_A[s],b_A[s],g_A[s]);
      }
    }
    if(inAlignment.EqualTo(" ")) {        // miniMC data
      printf("=== errors ==\n");
      Double_t eR0x=0,eR0y=0,eR0z=0,ea0=0,eb0=0,eg0=0;
      for(Int_t s=0;s<12;s++) 
       if(Sectors[s]) {
          dR0x[s]=rr0[0][s]-R0_A[0][s];
          dR0y[s]=rr0[1][s]-R0_A[1][s];
          dR0z[s]=rr0[2][s]-R0_A[2][s];
          da0[s]=ra[s]-a_A[s];
          db0[s]=rb[s]-b_A[s];
          dg0[s]=rg[s]-g_A[s];
          printf("sector=%2d",s);
          printf("  R0(%7.4f,%7.4f,%7.4f)  abg(%7.4f,%7.4f,%7.4f)\n",
                 dR0x[s],dR0y[s],dR0z[s],da0[s],db0[s],dg0[s]);
          eR0x+=dR0x[s]/ns,eR0y+=dR0y[s]/ns,eR0z+=dR0z[s]/ns,
          ea0+=da0[s]/ns,eb0+=db0[s]/ns,eg0+=dg0[s]/ns;
      }
      printf("=== errors shifts: (%7.4f,%7.4f,%7.4f  %7.4f,%7.4f,%7.4f)\n",
             eR0x,eR0y,eR0z,ea0,eb0,eg0);
      for(Int_t s=0;s<12;s++) {
        if(Sectors[s]) {
          printf("sector=%2d",s);
          printf("  R0(%7.4f,%7.4f,%7.4f)  abg(%7.4f,%7.4f,%7.4f)\n",
                 dR0x[s]-eR0x,dR0y[s]-eR0y,dR0z[s]-eR0z,
                 da0[s]-ea0,db0[s]-eb0,dg0[s]-eg0);
        }
      }
    }
    if(F>F0) {
      if(retC==0) 
      {
        printf("STOP! the iterative process does not converge. F=%f\n",F);
        cleanfiles(OutFile,AFile);
        return -1;
      }
      else {
        retC--;
        printf("<================ retN=%d ===========================================>f\n",retN);
       }
    }
    else {
      retC=retN;
    }
    if(eps<=Epsilon) {
      printf("====> epsilon limit is reached =====> eps=%f\n",eps);
      cleanfiles(OutFile,AFile);
      return 0;
    }
    if(loop==nLoop) {     //  last loop?
      printf("====> maximum loop number=%d is reached =====>\n",nLoop);
      cleanfiles(OutFile,AFile);
      return nLoop;
    }

    // find the next step to the minimum
//    geosec->PrintSectorAlignment(0);
///    miniTpcAlignment* al=new miniTpcAlignment();
    al->Init(*geosec, 0,0,AbsMinHits,Sectors,OutFile," ");
//    geosec->PrintSectorAlignment(0);
//    al->SetLaserMode(LM);
    al->SetPrintLevel(alignmentPL);
    al->alignment();
    al->results(dRx,dRy,dRz,da,db,dg);
    checkdiv(dRx,dRy,dRz,da,db,dg,Sectors);
///    delete al;
    F0=F;
  }
  printf("rm OutFile is %d    rm AFile is %d\n",remove(OutFile),remove(AFile));

  return -111;
}

void getAlignment(const char* afile) {
  TFile* file = new TFile(afile,"READ");
  TTree* algn = (TTree*)file->Get("alignment");
  algn->SetBranchAddress("R0_A",&R0_A);
  algn->SetBranchAddress("alpha_A",&a_A);
  algn->SetBranchAddress("beta_A",&b_A);
  algn->SetBranchAddress("gamma_A",&g_A);
  algn->GetEvent(0);
  file->Close();
  if(11==1){
  cout<<"===== getAlignment("<<afile<<") ======"<<endl;
  for(Int_t is=0;is<24;is++)
    printf("sector=%2d  R0(%11.4e,%11.4e,%11.4e)  abg(%11.4e,%11.4e,%11.4e)\n",
           is,R0_A[0][is],R0_A[1][is],R0_A[2][is],a_A[is],b_A[is],g_A[is]);
  cout<<"=================== =================="<<endl;
  }
}

void getRealAlignment(const char* afile) {
  TFile* file = new TFile(afile,"READ");
  TTree* algn = (TTree*)file->Get("real_alignment");
  algn->SetBranchAddress("R0_A",&R0_A);
  algn->SetBranchAddress("alpha_A",&a_A);
  algn->SetBranchAddress("beta_A",&b_A);
  algn->SetBranchAddress("gamma_A",&g_A);
  algn->GetEvent(0);
  delete algn;
  delete file;
//  file->Close();
}

void putAlignment(const char* afile) {
  TFile* file1 =  new TFile(afile,"RECREATE");
  TTree* algnw= new TTree("alignment","alignment_parameters");;
  algnw->Branch("R0_A",&R0_A,"R0_A[3][24]/D");
  algnw->Branch("alpha_A",&a_A,"a_A[24]/D");
  algnw->Branch("beta_A",&b_A,"b_A[24]/D");
  algnw->Branch("gamma_A",&g_A,"g_A[24]/D");
  algnw->Fill();
  algnw->Write();
  file1->Close();
  if(1==2){
  cout<<"===== putAlignment("<<afile<<") ======"<<endl;
  for(Int_t is=0;is<24;is++)
    printf("sector=%2d  R0(%11.4e,%11.4e,%11.4e)  abg(%11.4e,%11.4e,%11.4e)\n",
           is,R0_A[0][is],R0_A[1][is],R0_A[2][is],a_A[is],b_A[is],g_A[is]);
  cout<<"=================== =================="<<endl;
  }
}

Double_t getChi2(const char* infile) {
  Long64_t nhs;
  Double_t sum;
  TFile* file =  new TFile(infile,"READ");
  TTree* tchi2 = (TTree*)file->Get("Chi2");
  tchi2->SetBranchAddress("secNhits",&NhitSec0);
  tchi2->SetBranchAddress("secSumChi2",&SecChi2_0);
  tchi2->GetEvent(0);
  file->Close();
  sum=0;
  nhs=0;
  for(Int_t is=0;is<12;is++) {
    sum+=SecChi2_0[is];
    nhs+=NhitSec0[is];
    SecChi2_0[is]=(NhitSec0[is]>0)?SecChi2_0[is]/NhitSec0[is]:0;
  }
   Double_t Max=0;
   for(Int_t s=0;s<12;s++)
    if(SecChi2_0[s]>Max) {Max=SecChi2_0[s];maxS=s;}

  return (nhs>0)?sum/nhs:0;
}

Double_t getNhits(const char* infile) {
  Long64_t nhs;
  TFile* file =  new TFile(infile,"READ");
  TTree* tchi2 = (TTree*)file->Get("Chi2");
  tchi2->SetBranchAddress("secNhits",&NhitSec0);
  tchi2->SetBranchAddress("secSumChi2",&SecChi2_0);
  tchi2->GetEvent(0);
  file->Close();
  nhs=0;
  for(Int_t is=0;is<12;is++) {
    nhs+=NhitSec0[is];
  } 
  return nhs;
}

Double_t checkdiv(Double_t* dx,Double_t* dy,Double_t* dz,
                  Double_t* da,Double_t* db,Double_t* dg,Int_t* Sectors) {
  Double_t dX[3],dGX[3],dA[3],divX[3],divA[3];
    
  for(Int_t i=0;i<3;i++) {dGX[i]=0,dX[i]=0,dA[i]=0,divX[i]=0,divA[i]=0;}
  for(Int_t s=0;s<24;s++) {
    if(Sectors[s]) {
      dX[0]+=dx[s];dX[1]+=dy[s];dX[2]+=dz[s];
      divX[0]+=dx[s]*dx[s];divX[1]+=dy[s]*dy[s];divX[2]+=dz[s]*dz[s];
      dA[0]+=da[s];dA[1]+=db[s];dA[2]+=dg[s];
      divA[0]+=da[s]*da[s];divA[1]+=db[s]*db[s];divA[2]+=dg[s]*dg[s];
    }
  }
  for(Int_t i=0;i<3;i++) {
    dX[i]/=ns;dA[i]/=ns;divX[i]/=ns;divA[i]/=ns;
    divX[i]-=dX[i]*dX[i],
    divA[i]-=dA[i]*dA[i];
    divX[i]=TMath::Sqrt(divX[i]),
    divA[i]=TMath::Sqrt(divA[i]);
  }

  printf("dX  (%11.4e+-%11.4e   %11.4e+-%11.4e   %11.4e+-%11.4e)\n",
         dX[0],divX[0],dX[1],divX[1],dX[2],divX[2]);
  printf("dA  (%11.4e+-%11.4e   %11.4e+-%11.4e   %11.4e+-%11.4e)\n",
         dA[0],divA[0],dA[1],divA[1],dA[2],divA[2]);

  for(Int_t s=0;s<24;s++) {
    if(Sectors[s]) {
      dx[s]-=dX[0];dy[s]-=dX[1];dz[s]-=dX[2];
      for(Int_t k=0;k<3;k++) gdX[k]=dX[k];
    }
  }
  return 1;
}

void convert2miniMC(const char* InDataFile,const char* InAlignmentFile) {
};
