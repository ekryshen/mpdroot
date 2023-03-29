using ROOT::Math::cyl_bessel_i;
using TMath::ATan2;
using TMath::Cos;
using TMath::Pi;
using TMath::Sin;
using TMath::Sqrt;
using std::cout;
using std::cerr;
using std::endl;

double CalcResolution(double chi, double harmonic);
double GetChi(double res, double harmonic, int Niter);

void getResolution(std::string iFileName="pEP.root", std::string oFileName="pEP_resolution.root") 
{
   gSystem->Load("libMathMore.so");
   // Get input information
   TFile *fi = new TFile(iFileName.c_str(),"read");
   if (!fi) {
      cerr << "No input file was found!" << endl;
      return;
   }

   TProfile *mhCosFHCalNFHCalSAll = (TProfile*) fi->Get("mhCosFHCalNFHCalSAll");
   if (!mhCosFHCalNFHCalSAll) {
      cerr << "No histogram named mhCosFHCalNFHCalSAll was found in the ROOT file!" << endl;
      return;
   }
   TProfile *mhCosTPCNTPCSAll = (TProfile*) fi->Get("mhCosTPCNTPCSAll");
   if (!mhCosTPCNTPCSAll) {
      cerr << "No histogram named mhCosTPCNTPCSAll was found in the ROOT file!" << endl;
      return;
   }

   // Initialize output
   TFile *fo = new TFile(oFileName.c_str(),"recreate");

   // TProfile *mhRes1FHCalNFHCalSAll = (TProfile*) mhCosFHCalNFHCalSAll->Clone();
   // mhRes1FHCalNFHCalSAll->SetName("mhRes1FHCalNFHCalSAll");
   // TProfile *mhRes2FHCalNFHCalSAll = (TProfile*) mhCosFHCalNFHCalSAll->Clone();
   // mhRes2FHCalNFHCalSAll->SetName("mhRes2FHCalNFHCalSAll");
   // TProfile *mhRes1FHCalFullAll    = (TProfile*) mhCosFHCalNFHCalSAll->Clone();
   // mhRes1FHCalFullAll->SetName("mhRes1FHCalFullAll");
   // TProfile *mhRes2FHCalFullAll    = (TProfile*) mhCosFHCalNFHCalSAll->Clone();
   // mhRes2FHCalFullAll->SetName("mhRes2FHCalFullAll");
   // TProfile *mhRes2TPCNTPCSAll     = (TProfile*) mhCosTPCNTPCSAll->Clone();
   // mhRes2TPCNTPCSAll->SetName("mhRes2TPCNTPCSAll");

   TH1F *mhRes1FHCalNFHCalSAll = new TH1F("mhRes1FHCalNFHCalSAll","", 10, 0., 100.);
   TH1F *mhRes2FHCalNFHCalSAll = new TH1F("mhRes2FHCalNFHCalSAll","", 10, 0., 100.);
   TH1F *mhRes1FHCalFullAll    = new TH1F("mhRes1FHCalFullAll"   ,"", 10, 0., 100.);
   TH1F *mhRes2FHCalFullAll    = new TH1F("mhRes2FHCalFullAll"   ,"", 10, 0., 100.);
   TH1F *mhRes2TPCNTPCSAll     = new TH1F("mhRes2TPCNTPCSAll"    ,"", 10, 0., 100.);

   double res2, eres2, res, eres, chi, chiF;
   double resMin, resMax, chiMin, chiMax, chiMinF, chiMaxF, eresCalc, eresCalcF;
   // Calculate resolution for FHCal
   // Errors are estimated rather roughly from variations (mean+-error) 
   // of the original <cos(n(Psi1-Psi2))> correlation.
   // The relation between <cos(n(Psi1-Psi2))> and Rn is assumed to be unambiguous
   // and straightforward (i.e. with larger <cos(n(Psi1-Psi2))> we get larger Rn)
   for (int i=0; i<mhCosFHCalNFHCalSAll->GetNbinsX(); i++) {
      res2 = mhCosFHCalNFHCalSAll->GetBinContent(i+1);
      eres2= mhCosFHCalNFHCalSAll->GetBinError(i+1);
      eres = eres2/(2.*res2);
      res2 = fabs(res2);
      res  = Sqrt(res2);

      resMin = res - eres;
      resMax = res + eres;

      chi = GetChi(res, 1., 50);
      chiF= chi*Sqrt(2);

      if (resMin>0) {
         chiMin = GetChi(resMin, 1., 50);
         chiMinF= chiMin*Sqrt(2);
      }
      if (resMax<1) {
         chiMax = GetChi(resMax, 1., 50);
         chiMaxF= chiMax*Sqrt(2);
      }

      // Resolution for v1 from FHCal N/S
      res = CalcResolution(chi, 1.);
      mhRes1FHCalNFHCalSAll->SetBinContent(i+1, res);
      if (resMin>0) resMin = CalcResolution(chiMin, 1.);
      if (resMax<1) resMax = CalcResolution(chiMax, 1.);
      eres = (resMin>0 && resMax<1) ? fabs(resMax-resMin)/2. : 0.;
      mhRes1FHCalNFHCalSAll->SetBinError(i+1, eres);

      // Resolution for v2 from FHCal N/S
      res = CalcResolution(chi, 2.);
      mhRes2FHCalNFHCalSAll->SetBinContent(i+1, res);
      if (resMin>0) resMin = CalcResolution(chiMin, 1.);
      if (resMax<1) resMax = CalcResolution(chiMax, 1.);
      eres = (resMin>0 && resMax<1) ? fabs(resMax-resMin)/2. : 0.;
      mhRes2FHCalNFHCalSAll->SetBinError(i+1, eres);

      // Resolution for v1 from FHCal N+S
      res = CalcResolution(chiF, 1.);
      mhRes1FHCalFullAll->SetBinContent(i+1, res);
      if (resMin>0) resMin = CalcResolution(chiMin, 1.);
      if (resMax<1) resMax = CalcResolution(chiMax, 1.);
      eres = (resMin>0 && resMax<1) ? fabs(resMax-resMin)/2. : 0.;
      mhRes1FHCalFullAll->SetBinError(i+1, eres);

      // Resolution for v2 from FHCal N+S
      res = CalcResolution(chiF, 2.);
      mhRes2FHCalFullAll->SetBinContent(i+1, res);
      if (resMin>0) resMin = CalcResolution(chiMin, 1.);
      if (resMax<1) resMax = CalcResolution(chiMax, 1.);
      eres = (resMin>0 && resMax<1) ? fabs(resMax-resMin)/2. : 0.;
      mhRes2FHCalFullAll->SetBinError(i+1, eres);
   }

   // Calculate resolution for TPC N/S - eta-sub method
   for (int i=0; i<mhCosTPCNTPCSAll->GetNbinsX(); i++) {
      res2 = mhCosTPCNTPCSAll->GetBinContent(i+1);
      res2 = fabs(res2);
      res  = Sqrt(res2);
      eres2= mhCosTPCNTPCSAll->GetBinError(i+1);
      eres = eres2/(2.*res2);

      // Resolution for v2 from TPC N/S
      mhRes2TPCNTPCSAll->SetBinContent(i+1, res);
      mhRes2TPCNTPCSAll->SetBinError(i+1, eres);
   }

   // Add brief title to the output profiles
   mhRes1FHCalNFHCalSAll->SetTitle("R_{1} for FHCal N/S;Centrality, %;R_{1}");
   mhRes2FHCalNFHCalSAll->SetTitle("R_{2} for FHCal N/S;Centrality, %;R_{2}");
   mhRes1FHCalFullAll   ->SetTitle("R_{1} for FHCal F;Centrality, %;R_{1}");
   mhRes2FHCalFullAll   ->SetTitle("R_{2} for FHCal F;Centrality, %;R_{2}");
   mhRes2TPCNTPCSAll    ->SetTitle("R_{2} for TPC N/S;Centrality, %;R_{2}");

   fo->cd();
   mhRes1FHCalNFHCalSAll->Write();
   mhRes2FHCalNFHCalSAll->Write();
   mhRes1FHCalFullAll   ->Write();
   mhRes2FHCalFullAll   ->Write();
   mhRes2TPCNTPCSAll    ->Write();
   fo->Close();
}

double CalcResolution(double chi, double harmonic)
{
   Double_t con = Sqrt(Pi() / 2.) / 2.;
   Double_t arg = chi * chi / 4.;
   Double_t res =
      con * chi * exp(-arg) * (cyl_bessel_i(0.5 * (harmonic - 1), arg) + cyl_bessel_i(0.5 * (harmonic + 1), arg));
   return res;
}

double GetChi(double res, double harmonic, int Niter)
{
   Double_t chi   = 2.0;
   Double_t delta = 1.0;
   for (int i = 0; i < Niter; i++) {
      if (CalcResolution(chi, harmonic) < res) {
         chi = chi + delta;
      } else {
         chi = chi - delta;
      }
      delta = delta / 2.;
   }

   return chi;
}