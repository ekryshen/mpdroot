/**
 * \class MpdFemtoModelManager
 * \brief Manager for model studies
 *
 * The MpdFemtoModelManager class is a main helper class for femtoscopy
 * calculations. It manages weight generation, freeze-out coordinates
 * generation
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoModelManager_h
#define MpdFemtoModelManager_h

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseModelWeightGenerator.h"
#include "MpdFemtoBaseModelFreezeOutGenerator.h"
// Infrastructure
#include "MpdFemtoEnumeration.h"

//_________________
class MpdFemtoModelManager {
 public:
  /// Default constructor
  MpdFemtoModelManager();
  /// Copy constructor
  MpdFemtoModelManager(const MpdFemtoModelManager& manager);
  /// Assignment operator
  MpdFemtoModelManager& operator=(const MpdFemtoModelManager& manager);
  /// Destructor
  virtual ~MpdFemtoModelManager();

  /// Set pointer to the freeze-out coordinate generator
  void setFreezeOutGenerator(MpdFemtoBaseModelFreezeOutGenerator *foGen) {
    mFreezeOutGenerator = foGen;
  }
  /// Set pointer to the freeze-out coordinate generator
  void SetFreezeOutGenerator(MpdFemtoBaseModelFreezeOutGenerator *foGen) {
    setFreezeOutGenerator(foGen);
  }
  /// Set pointer to the femtoscopic weight generator
  void setWeightGenerator(MpdFemtoBaseModelWeightGenerator *weightGen) {
    mWeightGenerator = weightGen;
  }
  /// Set pointer to the femtoscopic weight generator
  void SetWeightGenerator(MpdFemtoBaseModelWeightGenerator *weightGen) {
    setWeightGenerator(weightGen);
  }
  /// Set status of the copy hidden info
  void createCopyHiddenInfo(bool aCopy = true);
  /// Set status of the copy hidden info
  void CreateCopyHiddenInfo(bool aCopy) {
    createCopyHiddenInfo(aCopy);
  }

  /// Return pointer to the freeze-out coordinate generator
  MpdFemtoBaseModelFreezeOutGenerator* freezeOutGenerator() {
    return mFreezeOutGenerator;
  }
  /// Return pointer to the freeze-out coordinate generator
  MpdFemtoBaseModelFreezeOutGenerator* FreezeOutGenerator() {
    return freezeOutGenerator();
  }
  /// Return pointer to the pair-weight generator
  MpdFemtoBaseModelWeightGenerator* weightGenerator() {
    return mWeightGenerator;
  }
  /// Return pointer to the pair-weight generator
  MpdFemtoBaseModelWeightGenerator* WeightGenerator() {
    return weightGenerator();
  }
  /// Return femtoscopic weight
  virtual double weight(MpdFemtoPair *aPair);
  /// Return femtoscopic weight
  virtual double Weight(MpdFemtoPair *aPair) {
    return weight(aPair);
  }

 protected:

  /// Pointer to freeze-out coordinate generator
  MpdFemtoBaseModelFreezeOutGenerator *mFreezeOutGenerator;
  /// Femtoscopic weight generator
  MpdFemtoBaseModelWeightGenerator *mWeightGenerator;
  /// Switch to turn on hidden-info generation
  bool mCreateCopyHiddenInfo;

 private:

  ClassDef(MpdFemtoModelManager, 1);
};

#endif // #define MpdFemtoModelManager_h
