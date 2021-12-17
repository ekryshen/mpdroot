/**********************************************************
** Track-based alignment using a Kalman filter technique **
** Kirill FOMENKO, LHE JINR, MPD Collaboration           **
***********************************************************/

#ifndef MPDTPCALIGNMENTKALMAN_H
#define MPDTPCALIGNMENTKALMAN_H

#include <random>
#include <cassert>

#include <FairTask.h>
#include <TClonesArray.h>
#include <TRotation.h>
#include <TH1.h>
#include <TFile.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

//#include "MpdTpcAlignmentParams.h"
#include "MpdTpcHit.h"
#include "MpdTpcSectorGeo.h"

#define DEBUG

class MpdTpcAlignmentKalman : public FairTask {
public:
   MpdTpcAlignmentKalman(Option_t *calib_tracks_type); // type = laser_mc/emulation/laser_data/cosmic_mu
   ~MpdTpcAlignmentKalman() override;

   // Task interfaces
   virtual InitStatus Init() override;
   virtual void       Exec(Option_t *) override;
   virtual void       Finish() override;
   virtual void       FinishEvent() override;
   void               SetPersistence(Bool_t val = false) { fPersistence = val; }

private:
   // dimensions
   static constexpr Int_t ndims = 3;             // coordinate space dim
   static constexpr Int_t npars = 3;             // alignment state vector dim
   static constexpr Int_t nipar = 2 * ndims;     // alignment vector begin row/column index in combined matrixes
   static constexpr Int_t nsums = nipar + npars; // combined state vector dim
   // static constexpr Int_t                                           nwmod = 3; // number of alignment pars need to be
   // corrected for weak mods

   // geometry TODO: params must be read from DB
   static constexpr Int_t  nSectors      = 24;
   static constexpr Int_t  nSectors_half = nSectors / 2;
   static constexpr double dPhi          = 2. * M_PI / nSectors_half;

   static constexpr double driftLength = 163.0;
   static constexpr double lowEdgeY    = 38.9;

   // track quality req.
   static constexpr Int_t nHitsPerSecMin      = 7;
   static constexpr Int_t nSectorsPerTrackMin = 2;

   // StatBuffer size
   static constexpr Int_t Awnd = 1000;

   // annealing factors
   static constexpr Int_t  nAnnSecTracks = 1000;
   static constexpr double ann0          = 1.E5;

   // cov matrix conditions
   static constexpr double Cmin = 5.E-5;
   static constexpr double Cini = 1.E-3;
   static constexpr double Qini = 1.E-4;

   // physics
   static constexpr double c_vel = 29.9792458; // cm/ns

   enum sums_names { s_x, s_y, s_z, s_vx, s_vy, s_vz, s_dx, s_dy, s_dz };
   enum pars_names { p_dx, p_dy, p_dz, p_dp, p_de, p_dk };

   struct Hit {
      double          time; //
      unsigned        sid;  // sector id
      Eigen::Vector3d raw;  // global coordinates, not corrected
      Eigen::Vector3d glb;  // global coordinates, corrected to the current alignment
      Eigen::Vector3d loc;  // local coordinates, not corrected
      Eigen::Vector3d fit;  // global coord. of the closest point on the fitted track
   };

   struct FittedTrack {
      FittedTrack() { reset(); }
      ~FittedTrack() { clear(); }

      void reset()
      {
         hits.clear();
         for (unsigned i = 0; i < nSectors; ++i) cov[i] = nullptr;
      }
      void clear()
      {
         hits.clear();
         for (unsigned i = 0; i < nSectors; ++i)
            if (cov[i]) {
               delete cov[i];
               cov[i] = nullptr;
            }
      }

      Int_t              side;          // -1 if track on the East side of the detector (z_glob < 0), 1 otherwise
      Int_t              nsec;          // number of sectors hit by the track's passage
      double             chi2;          // chi^2 of the track
      Eigen::Vector3d    org;           // starting point, global coordinates
      Eigen::Vector3d    dir;           // normalized dir vector at the starting point, global coordinates
      Eigen::Vector3d    org_mc;        // MC org
      Eigen::Vector3d    dir_mc;        // MC dir
      Eigen::Matrix3d   *cov[nSectors]; // covariances per sector (aka V)
      std::vector<Hit *> hits; // supposed to belong to the track (just refs to external hits, no need to take care)
   };

   // ring buffer to smooth params over size_max steps
   class StatBuffer {
   private:
      Eigen::Matrix<double, npars, 1> *data;
      size_t                           size_max, offset;
      bool                             data_ready;

   public:
      StatBuffer() : data(nullptr), size_max(0), offset(0), data_ready(false) {}
      ~StatBuffer()
      {
         if (data) delete[] data;
      }

      bool full_stats_available_() const { return data_ready; }

      StatBuffer &init(size_t size);
      StatBuffer &erase();
      StatBuffer &add(const Eigen::Matrix<double, npars, 1> &val);
      StatBuffer &stats_compute(Eigen::Matrix<double, npars, 1> &mean, Eigen::Matrix<double, npars, 1> &disp);
   };

   // converters
   Hit                                 *hit_(const MpdTpcHit *tpcHit) const;
   unsigned                             sid_(const Eigen::Vector3d &globXYZ) const;
   std::pair<Eigen::Vector3d, unsigned> glob2loc(const Eigen::Vector3d &globXYZ, const unsigned sid = nSectors) const;
   Eigen::Vector3d                      loc2glob(const Eigen::Vector3d &locXYZ, const unsigned &sid) const;

   // algo subroutines
   Eigen::Vector3d                      align(const unsigned &sid, const Eigen::Vector3d &p_loc,
                                              const Eigen::Matrix<double, npars, 1> *a = nullptr) const; // loc->glob, alignment applied
   std::pair<Eigen::Vector3d, unsigned> misalign(
      const Eigen::Vector3d                 &p_glb,
      const Eigen::Matrix<double, npars, 1> *a = nullptr) const; //  glob->loc, misalignment applied
   void matrixHInit(const unsigned &sid);

   bool sector_status_(const unsigned &sid);
   bool sector_state_init(const unsigned int &sid, const Eigen::Matrix<double, ndims, 1> &vel,
                          const Eigen::Matrix<double, ndims, 1> &pos, const FittedTrack *tr);
   void sector_state_revise(const unsigned &sid);
   void weak_modes_fix(Eigen::Matrix<double, npars, 1> *pars2fix);

   // workflow functions
   Int_t tracksFind();                                // returnes number of tracks found
   void  tracksProcess();                             // core function, alignment algo
   bool  trackParsCompute(FittedTrack &track_to_fit); // fitting procedure, track's params

   // per sector elements
   Eigen::Matrix<double, npars, 1>     s[nSectors]; // sector state vectors, to be found
   Eigen::Matrix<double, npars, npars> C[nSectors]; // cov(s), to be found

   Eigen::Matrix<double, nsums, nsums> A;           // propogation matrix
   Eigen::Matrix<double, nsums, nsums> Q;           // cov(A), includes scattering if needed
   Eigen::Matrix<double, ndims, nsums> H[nSectors]; // model matrixes
   Eigen::Matrix<double, ndims, ndims> V;           // cov(H), from coordinates' resolutions

   // global state matrixes
   Eigen::Matrix<double, nsums, 1>     sp; // corrected
   Eigen::Matrix<double, nsums, nsums> Cp; // corrected

   // service matrixes
   Eigen::Matrix<double, nsums, nsums> Ins;
   Eigen::Matrix<double, npars, npars> Inp;
   Eigen::Matrix<double, ndims, ndims> Ind;

   // data
   bool is_finished[nSectors];
   bool is_active[nSectors];

   Int_t events_processed;
   Int_t tracks_processed;

   Int_t  track_segments[nSectors];
   double ann[nSectors];
   double chi2[nSectors];

   StatBuffer pstats[nSectors];

   std::vector<Hit *>       hits;
   std::vector<FittedTrack> tracks;

   TClonesArray *fHits;
   char         *tracks_type;

#ifdef DEBUG
   Eigen::Matrix<double, npars, 1> s0[nSectors]; // initial alignment vectors
   Eigen::Matrix<double, npars, 1> st[nSectors]; // true alignment vectors (from MC)

   // MpdTpcAlignmentParams*                                           reader;
   MpdTpcSectorGeo *fSecGeo;
   TClonesArray    *fHitsAlignTrue;

   struct OutStruct {
      Double_t t[nSectors];
      Double_t ann[nSectors];
      Double_t nsg[nSectors];
      Double_t nht[nSectors];
      Double_t dx[nSectors];
      Double_t dy[nSectors];
      Double_t dz[nSectors];
      Double_t dp[nSectors];
      Double_t de[nSectors];
      Double_t dk[nSectors];
      Double_t Cx[nSectors];
      Double_t Cy[nSectors];
      Double_t Cz[nSectors];
      Double_t Cp[nSectors];
      Double_t Ce[nSectors];
      Double_t Ck[nSectors];
      // Double_t Vx[nSectors];
      // Double_t Vy[nSectors];
      // Double_t Vz[nSectors];
   } ost;

   TFile *out_file;
   TTree *out_tree;

   TH1I *hits_per_trk;
   TH1I *hits_per_seg;
   TH1I *segs_per_trk;
   TH1I *hits_per_sec[nSectors];
   TH1I *segs_per_sec[nSectors];

   std::mt19937 *rand_gen;

   double rnd_flat(const double min, const double max) const; // NB: [min, max)
   double rnd_norm(const double mean, const double sigma) const;
   double rnd_exp(const double mean) const;

   unsigned printout_period;

   void misalignment_commit(unsigned nHits);
   void event_log_flush();
   void matrix_print(unsigned sid);

   // void                                                             hit2vector(const Hit* hit, TVector3& globXYZ,
   // TVector3& locXYZ) const; std::pair<Eigen::Vector3d, unsigned>                             glob2loc(const TVector3&
   // globXYZ) const; Eigen::Vector3d                                                  loc2glob(const TVector3& locXYZ,
   // const unsigned& sid) const; TVector3                                                         toTV3(const
   // Eigen::Vector3d& v) const { return TVector3(v(0,0), v(1,0), v(2,0)); }

// public:
// void                                                             SetAlignParamsFile(std::string path) { if(reader)
// reader->SetAlignParamsFile(path); }
#endif

   Bool_t  fPersistence;
   clock_t exec_time;

   ClassDefOverride(MpdTpcAlignmentKalman, 0)
};

#endif
