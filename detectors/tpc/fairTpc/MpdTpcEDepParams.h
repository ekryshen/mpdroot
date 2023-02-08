#ifndef MPDTPCEDEPPARAMS_H
#define MPDTPCEDEPPARAMS_H

/// \ingroup sim
/// \class MpdTpcEDepParams
/// \brief Parameters for dE/dx simulation in MPD TPC
///
/// \author Igor Rufanov, LHEP JINR Dubna - development
/// \author Alexander Zinchenko, LHEP JINR Dubna - porting to MpdRoot
/// 18.11.2019

#include <TObject.h>
//#include <TVector3.h>

#include <vector>

class MpdTpcEDepParams : public TObject {

public:
   static MpdTpcEDepParams *Instance(); ///< get singleton instance

   Double_t GetCollisionDensity(Double_t log10bg);
   Double_t GetRandEnergy(Double_t log10bg);
   Double_t GetMinLimit() const { return fCollDensL10BgMin; }
   Double_t GetMaxLimit() const { return fCollDensL10BgMax; }
   Double_t GetEloss(Double_t log10bg, Double_t charge, Double_t step);

protected:
   MpdTpcEDepParams() { Init(); }    ///< Default ctor
   virtual ~MpdTpcEDepParams() { ; } ///< Destructor

private:
   void Init();
   // Automatic deletion when application exits
   static void DestroyInstance()
   {
      if (fgTpcEdepParams) delete fgTpcEdepParams;
   }

   static MpdTpcEDepParams *fgTpcEdepParams; //! pointer to Singleton instance

   // static const Int_t fgkNsect = 12; // number of TPC sectors
   Int_t     fNCollDensBgBins;   // number of steps in log10(bet*gam)
   Double_t *fCollDens;          // array of cluster/cm values
   Double_t  fCollDensL10BgMin;  // min log10(bet*gam)
   Double_t  fCollDensL10BgMax;  // max log10(bet*gam)
   Double_t  fCollDensL10BgStep; // step in log10(bet*gam)

   Int_t fNEneBgBins;   // number of steps in log10(bet*gam)
   Int_t fNEneProbBins; // number of steps in probability
   // Double_t fEne[ TpcEDepParams::fNEneBgBins][TpcEDepParams::fNEneProbBins];
   std::vector<std::vector<Double_t>> fEne;          // energy loss vs prob vs log10(bet*gam)
   Double_t                           fEneL10BgMin;  // min log10(bet*gam)
   Double_t                           fEneL10BgMax;  // max log10(bet*gam)
   Double_t                           fEneL10BgStep; // step in log10(bet*gam)
   Double_t                           fEneProbStep;  // step in probability

   ClassDef(MpdTpcEDepParams, 0);
};

//__________________________________________________________________________

#endif
