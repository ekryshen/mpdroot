/**
 * \class MpdFemtoModelGausLCMSFreezeOutGenerator
 * \brief Generates freeze-out coordinates using Gaussian profile
 *
 * Freeze-out coordinates generator, generating a 3D gaussian ellipsoid in LCMS
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoModelGausLCMSFreezeOutGenerator_h
#define MpdFemtoModelGausLCMSFreezeOutGenerator_h

// MpdFemtoMaker headers
#include "MpdFemtoBaseModelFreezeOutGenerator.h"

// ROOT headers
#include "TRandom3.h"

//_________________
class MpdFemtoModelGausLCMSFreezeOutGenerator : public MpdFemtoBaseModelFreezeOutGenerator {
 public:
  /// Default constructor
  MpdFemtoModelGausLCMSFreezeOutGenerator();
  /// Copy construcotr
  MpdFemtoModelGausLCMSFreezeOutGenerator(const MpdFemtoModelGausLCMSFreezeOutGenerator &aModel);
  /// Assignment operator
  MpdFemtoModelGausLCMSFreezeOutGenerator& operator=(const MpdFemtoModelGausLCMSFreezeOutGenerator &aModel);
  /// Default destructor
  virtual ~MpdFemtoModelGausLCMSFreezeOutGenerator();

  /// Generate freeze-out coordinate (x,y,z,t)
  virtual void generateFreezeOut(MpdFemtoPair *aPair);
  /// Generate freeze-out coordinate (x,y,z,t)
  virtual void GenerateFreezeOut(MpdFemtoPair *aPair) {
    generateFreezeOut(aPair);
  }

  /// Set outward component
  void setSizeOut(const double& aSizeOut) {
    mSizeOut = aSizeOut;
  }
  /// Set sideward component
  void setSizeSide(const double& aSizeSide) {
    mSizeSide = aSizeSide;
  }
  /// Set longitudinal component
  void setSizeLong(const double& aSizeLong) {
    mSizeLong = aSizeLong;
  }

  /// Return outward component
  double sizeOut() const {
    return mSizeOut;
  }
  /// Return sideward component
  double sizeSide() const {
    return mSizeSide;
  }
  /// Return longitudinal component
  double sizeLong() const {
    return mSizeLong;
  }

  /// Clone generator
  virtual MpdFemtoBaseModelFreezeOutGenerator* clone() const {
    return generator();
  }

 protected:
  /// Size of the source in the out direction
  double mSizeOut;
  /// Size of the source in the side direction
  double mSizeSide;
  /// Size of the source in the long direction
  double mSizeLong;

 private:

  MpdFemtoBaseModelFreezeOutGenerator* generator() const {
    MpdFemtoBaseModelFreezeOutGenerator* tModel = new MpdFemtoModelGausLCMSFreezeOutGenerator(*this);
    return tModel;
  }

  ClassDef(MpdFemtoModelGausLCMSFreezeOutGenerator, 1);
};

#endif // MpdFemtoModelGausLCMSFreezeOutGenerator_h
