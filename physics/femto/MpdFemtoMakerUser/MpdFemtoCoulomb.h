/**
 * \class MpdFemtoCoulomb
 * \brief Class that performs Coulomb correction
 *
 *  This is a Coulomb correction class which
 *  1. Reads in the dat from a file
 *  2. Performs a linear interpolation in R and creates any array of interpolations
 *  3. Interpolates in eta and returns the Coulomb correction to user
 */

#ifndef MpdFemtoCoulomb_h
#define MpdFemtoCoulomb_h

// MpdFemtoMaker headers
#include "MpdFemtoPair.h"

// ROOT headers
#include "TH1D.h"
#include "TH3D.h"

//_________________

class MpdFemtoCoulomb {
public:
   /// Default constructor
   MpdFemtoCoulomb();
   /// Constructor with parameters
   MpdFemtoCoulomb(const char *readFile, const double &radius, const double &charge);
   /// Copy constructor
   MpdFemtoCoulomb(const MpdFemtoCoulomb &copy);
   /// Assignement operator
   MpdFemtoCoulomb &operator=(const MpdFemtoCoulomb &copy);
   /// Destructor
   virtual ~MpdFemtoCoulomb();

   /// Set radius for calculations
   void setRadius(const double &radius);
   /// Set output file name
   void setFile(const char *readFile);
   /// Set charge product
   void setChargeProduct(const double &charge);

   /// Radius

   double radius() { return mRadius; }
   /// These have different names so eta/Qinv don't confuse the compiler
   double coulombCorrect(const double &eta);
   /// These have different names so eta/Qinv don't confuse the compiler
   double coulombCorrect(const double &eta, const double &radius);
   /// Return Coulomb corrected pair

   double coulombCorrect(const MpdFemtoPair *pair) { return coulombCorrect(eta(pair)); }
   /// Return Coulomb corrected pair

   double coulombCorrect(const MpdFemtoPair *pair, const double &radius) { return coulombCorrect(eta(pair), radius); }
   /// Return Coulomb corrected value
   double coulombCorrect(const double &mass, const double &charge, const double &radius, const double &qInv);
   /// Return histogram with the Coulomb correction
   TH1D *correctionHistogram(const double &mass1, const double &mass2, const int &nBins, const double &low,
                             const double &high);
   TH1D *correctionHistogram(const TH1D *, const double);
   TH3D *correctionHistogram(const TH3D *, const double);

private:
   /// Calculates eta
   double eta(const MpdFemtoPair *pair);
   /// Creates look-up table
   void createLookupTable(const double &radius);
   /// File to interpolate corrections from
   const char *mFile;
   /// Radius from previous iteration
   double mRadius;
   /// Charge product of particles
   double mZ1Z2;
   /// Interpolated Coulomb correction table
   double mEta[1000];
   /// Interpolated Coulomb correction table
   double mCoulomb[1000];
   /// Number of Eta's in lookup-table
   int mNLines;

   ClassDef(MpdFemtoCoulomb, 0);
};

#endif
