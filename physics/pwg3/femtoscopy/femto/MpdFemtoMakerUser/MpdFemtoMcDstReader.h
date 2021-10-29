/**
 * \class MpdFemtoMcDstReader
 * \brief Reader for the McDst format that takes McDstReader
 *
 * The reader class for the McDst format. It reads McDst and
 * converts data to the internal MpdFemtoMaker structure (MpdFemtoEvent,
 * MpdFemtoTrack, etc).
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoMcDstReader_h
#define MpdFemtoMcDstReader_h

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseEventReader.h"
// Infrastructure
#include "MpdFemtoEvent.h"
#include "MpdFemtoTrack.h"
#include "MpdFemtoV0.h"
#include "MpdFemtoKink.h"
#include "MpdFemtoXi.h"
#include "MpdFemtoEnumeration.h"
#include "MpdFemtoString.h"

// ROOT headers
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"
#include "TRandom3.h"
#include "TMatrixTLazy.h"

// McDst headers
#include "MpdMcDst.h"
#include "MpdMcEvent.h"
#include "MpdMcParticle.h"
#include "MpdMcArrays.h"
#include "MpdMcDstReader.h"

//_________________
class MpdFemtoMcDstReader : public MpdFemtoBaseEventReader {
 public:
  /// Default constructor
  MpdFemtoMcDstReader();
  /// Parametrized constructor
  /// \param mcDstReader Takes pointer to McDstReader
  /// \param debug Set debut value (see the MpdFemtoBaseEventReader class)
  MpdFemtoMcDstReader(McDstReader* mcDstReader, int debug = 1);
  /// Copy constructor
  MpdFemtoMcDstReader(const MpdFemtoMcDstReader& copy);
  /// Assignment operator
  MpdFemtoMcDstReader& operator=(const MpdFemtoMcDstReader& copy);
  /// Destructor
  virtual ~MpdFemtoMcDstReader();

  /// Prepare report
  virtual MpdFemtoString report();
  /// Return instance of MpdFemtoEvent with filled information
  virtual MpdFemtoEvent* returnHbtEvent();

  /// Set McDstReader
  void setMcDstReader(McDstReader* reader) {
    mMcDstReader = reader;
  }
  
  /// Set magnetic field value in kilogauss (default is 0.5T)
  void setMagneticField(const float& field) {
    mMagField = field;
  }
  
  /// Set magnetic field value in kilogauss (default is 0.5T)
  void SetMagneticField(const float& field) {
    setMagneticField(field);
  }

  /// Set event plane rotation
  void setRotateEventPlane(const bool& rotate) {
    mDoRotate = rotate;
  }
  
  /// Set event plane rotation
  void SetRotateEventPlane(const bool& rotate) {
    setRotateEventPlane(rotate);
  }
  
  /// Set event plane angle range [min, max]
  void setEventPlaneAngleRange(const float& low, const float& hi) {
    mPhi[0] = low;
    mPhi[1] = hi;
  }
  
  /// Set event plane angle range [min, max]
  void SetEventPlaneAngleRange(const float& low, const float& hi) {
    setEventPlaneAngleRange(low, hi);
  }
  
  /// Set event plane resolution
  void setEventPlaneResolution(const float& res) {
    mEventPlaneResolution = res;
  }

  /// Set event plane resolution
  void SetEventPlaneResolution(const float& res) {
    setEventPlaneResolution(res);
  }

 private:

  /// Simple dE/dx estimation (not for physics analysis)
  double dEdxMean(Double_t mass, Double_t momentum);

  /// Pointer to the McDstReader
  McDstReader *mMcDstReader;

  /// Magnetic field in kilogauss (default is set to 0.5 T)
  float mMagField;

  /// Current event number
  //int mCurrentEvent;

  /// Pointer MpdFemtoEvent that should be return to MpdFemtoAnalysis
  MpdFemtoEvent *mHbtEvent;

  /// Use event plane rotation
  bool mDoRotate;
  /// Randomizer
  TRandom3 *mRandom;
  /// Phi range to generate event plane angle ( default: [-pi;pi] )
  float mPhi[2];
  /// Event plane resolution (default: ideal resolution, i.e. unity)
  float mEventPlaneResolution;

  /// Reference multiplicity (STAR: charged particles with pT>0.15, |eta|<0.5)
  unsigned short mRefMult;
  /// Reference multiplicity of positive particles (STAR: charged particles with pT>0.15, |eta|<0.5)
  unsigned short mRefMultPos;
  /// Reference multiplicity (STAR: charged particles with pT>0.15, |eta|<1)
  unsigned short mRefMult2;
  /// Reference multiplicity of positive particles (STAR: charged particles with pT>0.15, |eta|<1)
  unsigned short mRefMult2Pos;
  /// Sphericity esitmated by particles with  |eta|<0.5, pT>0.15
  float mSphericity;
  /// Sphericity esitmated by particles with  |eta|<1, pT>0.15
  float mSphericity2;

  /// Number of events passed
  long int mEventsPassed;

  /// Matrix for sphericity estimation (pT>0.15 && |eta|<0.5)
  TMatrixTSym<double> *mMatrix;
  /// Matrix for sphericity2 estimation (pT>0.15 && |eta|<1.0)
  TMatrixTSym<double> *mMatrix2;


  ClassDef(MpdFemtoMcDstReader, 3)
};

#endif // MpdFemtoMcDstReader_h
