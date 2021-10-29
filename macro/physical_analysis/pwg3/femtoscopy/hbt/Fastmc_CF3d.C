#include <TH3.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TString.h>
#include <iostream>
#include <TMath.h>


//-----------only FASTMC ----(HBT)----------

//3D -Correlation Functions for pions
//You have to compile it with using MakeFile
// make

// intervals of transverse momentum pt
 const int NPT = 3 ;
// TString binnames [NPT] = {"01-02","02-03","03-04","04-05","055-075","075-095"} ;
 TString binnames [NPT] = {"0-01","01-02","02-04"} ;
// double pt_low [NPT] = {0.1, 0.2, 0.3, 0.4, 0.55, 0.75} ;
 double pt_low [NPT] = {0., 0.1, 0.2} ;
// double pt_high [NPT] = {0.2, 0.3, 0.4, 0.5, 0.75, 0.95} ;
 double pt_high [NPT] = {0.1, 0.2, 0.4} ;
 double dpl = 0.5 ;
 const double QMAX = 0.08 ;
//const double QMAX = 0.1 ;

 const int NEVPACK = 50, NBINS = 8 ;

int main()
{

// FASTMC output file :
  TFile *f = new TFile("RunOutput.root") ;
// new ROOT-file:
  TFile outputFile("CFs.root", "RECREATE");

   TH3F *hCFnom [NPT];
   TH3F *hCFden [NPT];
   TH3F *hCF [NPT];
   
   TH1D *hCFnomL[NPT];
   TH1D *hCFdenL[NPT];
   TH1D *hCFL[NPT];

   TH1D *hCFnomS[NPT];
   TH1D *hCFdenS[NPT];
   TH1D *hCFS[NPT];    
   
   TH1D *hCFnomO[NPT];
   TH1D *hCFdenO[NPT];
   TH1D *hCFO[NPT];
   
   TH1D *hCFnomD[NPT];
   TH1D *hCFdenD[NPT];
   TH1D *hCFD[NPT];
   
   TH1D *hCFW[NPT];
    
   
   for(int ipt=0; ipt<NPT; ipt++){

     hCFnom [ipt] = new TH3F("hCFnom"+binnames[ipt],"hCFnom"+binnames[ipt],NBINS, 0., QMAX,NBINS, 0., QMAX,NBINS, 0., QMAX);
     hCFden [ipt] = new TH3F("hCFden"+binnames[ipt],"hCFden"+binnames[ipt],NBINS, 0., QMAX,NBINS, 0., QMAX,NBINS, 0., QMAX);
     hCF [ipt] = new TH3F("hCF"+binnames[ipt],"hCF"+binnames[ipt],NBINS, 0., QMAX,NBINS, 0., QMAX,NBINS, 0., QMAX);
     
       hCFnomL [ipt]  = new TH1D("hCFnomL"+binnames[ipt],"hCFnomL"+binnames[ipt],NBINS,0.,QMAX);
       hCFdenL [ipt]  = new TH1D("hCFdenL"+binnames[ipt],"hCFdenL"+binnames[ipt],NBINS,0.,QMAX);
       hCFL [ipt]  = new TH1D("hCFL"+binnames[ipt],"hCFL"+binnames[ipt],NBINS,0.,QMAX);  
       hCFnomS [ipt]  = new TH1D("hCFnomS"+binnames[ipt],"hCFnomS"+binnames[ipt],NBINS,0.,QMAX);
       hCFdenS [ipt]  = new TH1D("hCFdenS"+binnames[ipt],"hCFdenS"+binnames[ipt],NBINS,0.,QMAX);
       hCFS [ipt]  = new TH1D("hCFS"+binnames[ipt],"hCFS"+binnames[ipt],NBINS,0.,QMAX);  
       hCFnomO [ipt]  = new TH1D("hCFnomO"+binnames[ipt],"hCFnomO"+binnames[ipt],NBINS,0.,QMAX);
       hCFdenO [ipt]  = new TH1D("hCFdenO"+binnames[ipt],"hCFdenO"+binnames[ipt],NBINS,0.,QMAX);
       hCFO [ipt]  = new TH1D("hCFO"+binnames[ipt],"hCFO"+binnames[ipt],NBINS,0.,QMAX);     
       hCFnomD [ipt]  = new TH1D("hCFnomD"+binnames[ipt],"hCFnomD"+binnames[ipt],NBINS,0.,QMAX);
       hCFdenD [ipt]  = new TH1D("hCFdenLD"+binnames[ipt],"hCFdenD"+binnames[ipt],NBINS,0.,QMAX);
       hCFD [ipt]  = new TH1D("hCFD"+binnames[ipt],"hCFD"+binnames[ipt],NBINS,0.,QMAX); 
        hCFW [ipt]  = new TH1D("hCFW"+binnames[ipt],"hCFW"+binnames[ipt],NBINS,0.,QMAX);              
	   }

  const Int_t kMax = 1000000;

  Int_t npart[NEVPACK] ;
  Int_t  *pdg = new Int_t [kMax];
  Int_t   *Mpdg = new Int_t [kMax];
  Float_t *Px = new Float_t [kMax];
  Float_t *Py = new Float_t [kMax];
  Float_t *Pz = new Float_t [kMax];
  Float_t *E = new Float_t [kMax];
  Float_t *X = new Float_t [kMax];
  Float_t *Y = new Float_t [kMax];
  Float_t *Z = new Float_t [kMax];
  Float_t *T = new Float_t [kMax];

  Info("CF3d","arrays created\n");

  TTree *td = (TTree*)f->Get("td");
  Int_t nevents = td->GetEntries();
  Info("CF3d.C", "Nevents %d, ", nevents);

 TLorentzVector m1 ;
 TLorentzVector m2 ;
 TVector3 vlong ;

  for (Int_t k=0;k<nevents;k+=NEVPACK)
  {
        std::cout<<"Test k: "<<k<<std::endl;
  if (k%10 == 0) Info("","[%i]\n",k);

  int pcount = 0 ;

// take a group of NEVPACK particles
  for(int kk=0; kk<NEVPACK;kk++)
  {
    td->SetBranchAddress("npart",npart+kk);
    td->SetBranchAddress("Px",Px+pcount);
    td->SetBranchAddress("Py",Py+pcount);
    td->SetBranchAddress("Pz",Pz+pcount);
    td->SetBranchAddress("E",E+pcount);
    td->SetBranchAddress("X",X+pcount);
    td->SetBranchAddress("Y",Y+pcount);
    td->SetBranchAddress("Z",Z+pcount);
    td->SetBranchAddress("T",T+pcount);
    td->SetBranchAddress("pdg",pdg+pcount);
    td->SetBranchAddress("Mpdg",Mpdg+pcount);

    td->GetEntry(k+kk);

    pcount += npart[kk] ;
  }

     for (Int_t i=0;i<pcount;i++)
       if(pdg[i]==211){
	 // std::cout<<"Test i: "<<i<<std::endl;

        for (Int_t j=i+1; j<pcount; j++)
          if(pdg[j]==211){

            Float_t eta1 = -0.5*TMath::Log((E[i]+Pz[i])/(E[i]-Pz[i]+0.00001));
            Float_t eta2 = -0.5*TMath::Log((E[j]+Pz[j])/(E[j]-Pz[j]+0.00001));

	    // if (TMath::Abs(Pz[i])<dpl && TMath::Abs(Pz[j])<dpl){
	    if(abs(eta1)<0.5 && abs(eta2)<0.5){
		{
		m1.SetXYZT(Px[i],Py[i],Pz[i],E[i]);
		m2.SetXYZT(Px[j],Py[j],Pz[j],E[j]);
		vlong.SetXYZ(0., 0., (Pz[i] + Pz[j]) / (E[i] + E[j]));
		m1.Boost(-vlong);
		m2.Boost(-vlong);
              Float_t kt =  0.5 * TMath::Hypot(m1.X() + m2.X(), m1.Y() + m2.Y());
               for(int ipt=0; ipt<NPT; ipt++) // pt loop
                 if (kt>pt_low[ipt] && kt<pt_high[ipt]) {

               Float_t pXsum = m1.X() + m2.X();
               Float_t pYsum = m1.Y() + m2.Y();
               Float_t pXdif = m1.X() - m2.X();
               Float_t pYdif = m1.Y() - m2.Y();
               Float_t qOut = 0.5 * (pXsum * pXdif + pYsum * pYdif) / kt;

               Float_t qSide = (m1.X() * m2.Y() - m1.Y() * m2.X()) / kt;

               Float_t qLong = m1.Z() - m2.Z() ;
// ----------Test:
              Float_t dd =  1.*TMath::Sqrt((qOut*qOut)+(qSide*qSide)+(qLong*qLong));
	      
		 Float_t wCos = 1 + TMath::Cos(5.068423*((Px[i]-Px[j])*(X[i]-X[j])
                  + (Py[i]-Py[j])*(Y[i]-Y[j]) + (Pz[i]-Pz[j])*(Z[i]-Z[j])
                  - (T[i]-T[j])*(E[i]-E[j])));
                       //for pions
                        if (qOut > 0. && qOut < QMAX && qSide > 0. && qSide <QMAX && qLong > 0. && qLong <QMAX) {
                            hCFnom[ipt]->Fill(qOut,qSide,qLong, wCos);
                            hCFden[ipt]->Fill(qOut,qSide,qLong, 1.);
			    
			     hCFnomL[ipt]->Fill(qLong, wCos);
                             hCFdenL[ipt]->Fill(qLong, 1.);
			     hCFnomS[ipt]->Fill(qSide, wCos);
                             hCFdenS[ipt]->Fill(qSide, 1.);			     
			     hCFnomO[ipt]->Fill(qOut, wCos);
                             hCFdenO[ipt]->Fill(qOut, 1.);
			     hCFnomD[ipt]->Fill(dd, wCos);
                             hCFdenD[ipt]->Fill(dd, 1.);	
                             hCFW[ipt]->Fill(wCos, 1.);			     		    
                          }
                    } // pt loop
		}
                } //if the same particles

         } //for part2

     } //for part1

  } //for event


	for(int ipt=0; ipt<NPT; ipt++)
		{
   hCFnom[ipt]->Sumw2();
   hCFden[ipt]->Sumw2();

    hCFnomL[ipt]->Sumw2();
    hCFdenL[ipt]->Sumw2();
    hCFnomS[ipt]->Sumw2();
    hCFdenS[ipt]->Sumw2();    
     hCFnomO[ipt]->Sumw2();
    hCFdenO[ipt]->Sumw2();
    hCFnomD[ipt]->Sumw2();
    hCFdenD[ipt]->Sumw2();   
   hCFW[ipt]->Sumw2();     
    
    
   
   hCF[ipt]->Divide(hCFnom[ipt],hCFden[ipt], 1., 1.);
   hCF[ipt]->Write();
   
    hCFL[ipt]->Divide(hCFnomL[ipt],hCFdenL[ipt], 1., 1.);
    hCFL[ipt]->Write();
     hCFS[ipt]->Divide(hCFnomS[ipt],hCFdenS[ipt], 1., 1.);
    hCFS[ipt]->Write();   
    hCFO[ipt]->Divide(hCFnomO[ipt],hCFdenO[ipt], 1., 1.);
    hCFO[ipt]->Write();
     hCFD[ipt]->Divide(hCFnomD[ipt],hCFdenD[ipt], 1., 1.);
    hCFD[ipt]->Write(); 
     hCFW[ipt]->Write();   
    
   delete hCF[ipt] ;
   delete hCFnom[ipt] ;
   delete hCFden[ipt] ;
     delete hCFL[ipt] ;
   delete hCFnomL[ipt] ;
   delete hCFdenL[ipt] ; 
     delete hCFS[ipt] ;
   delete hCFnomS[ipt] ;
   delete hCFdenS[ipt] ;
      delete hCFO[ipt] ;
   delete hCFnomO[ipt] ;
   delete hCFdenO[ipt] ; 
     delete hCFD[ipt] ;
   delete hCFnomD[ipt] ;
   delete hCFdenD[ipt] ;  
         delete hCFW[ipt] ;
		}
   delete td ;

 outputFile.Write();


}
