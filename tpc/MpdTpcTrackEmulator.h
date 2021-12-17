#ifndef MPDTPCTRACKEMULATOR_H
#define MPDTPCTRACKEMULATOR_H

#include <random>
#include <Eigen/Core>

#include <FairTask.h>
#include <TClonesArray.h>

#include "MpdTpcHit.h"
#include "MpdTpcSectorGeo.h"
//#include "MpdTpcAlignmentParams.h"

class MpdTpcTrackEmulator : public FairTask {
public:
   // TODO: must be read from DB
   static constexpr double L_max  = 163.0; // drift length
   static constexpr double R_max  = 123.;  // outer radius
   static constexpr double R_min  = 38.9;  // pad planes low edge Y
   static constexpr double R2_max = R_max * R_max;
   static constexpr double R2_min = R_min * R_min;
   static constexpr double dPhi   = M_PI / 6.;

   MpdTpcTrackEmulator(Option_t *emu_tracks_type);
   virtual ~MpdTpcTrackEmulator() override { delete tracks_type; }

   virtual InitStatus Init() override;
   virtual void       Exec(Option_t *opt) override;
   virtual void       Finish() override
   {
      if (rand_gen) delete rand_gen;
   }

   MpdTpcTrackEmulator &stddev_set(const Double_t val)
   {
      stddev = val;
      return *this;
   }
   MpdTpcTrackEmulator &mean_free_path_set(const Double_t val)
   {
      mfpath = val;
      return *this;
   }
   MpdTpcTrackEmulator &hits_resolution_set(const Double_t val)
   {
      hitres = val;
      return *this;
   }
   MpdTpcTrackEmulator &velocity_set(const Double_t val)
   {
      velocity = val * 1.E-9;
      return *this;
   } // cm/s -> cm/ns
   MpdTpcTrackEmulator &tracks_per_event_set(const Int_t val)
   {
      ntrk_per_ev = val;
      return *this;
   }

private:
   void     vertexes_emulate();
   unsigned track_pass(Eigen::Vector3d const &org, Eigen::Vector3d const &dir);

   double rnd_flat(const double min, const double max) const; // NB: [min, max)
   double rnd_norm(const double mean, const double sigma) const;
   double rnd_exp(const double mean) const;

   char                                      *tracks_type;
   std::vector<std::pair<TVector3, Double_t>> hits;

   Double_t stddev;      // hit-to-track standart deviation. cm } have to
   Double_t mfpath;      // mean free path, cm                  } be set
   Double_t hitres;      // hits resolution, cm                 } explicitly
   Double_t velocity;    // velocity, cm/ns
   Int_t    ntrk_per_ev; // number of tracks in each event

   std::mt19937 *rand_gen;

   TClonesArray    *fHits;
   MpdTpcSectorGeo *fSecGeo;

   ClassDefOverride(MpdTpcTrackEmulator, 0)
};

#endif
